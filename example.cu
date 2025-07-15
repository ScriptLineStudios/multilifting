#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include "cubiomes/rng.h"
#include "cubiomes/finders.h"

typedef enum {
    STRUCTURE, BIOME, PILLAR
} ConstraintType;

typedef struct {
    ConstraintType tp; 

    int type;
    int x, z;
    
    int chunk_range, offset_x, offset_z, reg_x, reg_z;
    uint64_t salt;
} Constraint;

#define STRUCTURE_CONSTRAINT(t, _x, _z) (Constraint){.tp=STRUCTURE, .type=t, .x=_x, .z=_z}
#define BIOME_CONSTRAINT(t, _x, _z) (Constraint){.tp=BIOME, .type=t, .x=_x, .z=_z}
#define PILLAR_CONSTRAINT(seed) (Constraint){.tp=PILLAR, .type=seed}

typedef struct {
    Constraint *cons;
    size_t length;
} Constraints;

void parse_constraints(Constraints constraints) {
    if (constraints.length < 1) {
        fprintf(stderr, "ERROR: You need at least one constraint!");
        exit(1);
    }

    StructureConfig conf;
    int total_bits = 0;
    for (size_t i = 0; i < constraints.length; i++) {
        Constraint *c = &constraints.cons[i];
        if (c->tp != STRUCTURE) { // only structure constraints need to be parsed :D
            continue;
        }
 
        getStructureConfig(c->type, MC_1_16, &conf);
        c->salt = conf.salt;
        c->chunk_range = conf.chunkRange;
        
        int cx = c->x >> 4;
        int cz = c->z >> 4;

        int _cx = cx < 0 ? cx - conf.regionSize + 1 : cx;
        int _cz = cz < 0 ? cz - conf.regionSize + 1 : cz;

        c->reg_x = floor(_cx / conf.regionSize);
        c->reg_z = floor(_cz / conf.regionSize);

        c->offset_x = cx - (c->reg_x * conf.regionSize);
        c->offset_z = cz - (c->reg_z * conf.regionSize);

        total_bits += (int)log2((conf.chunkRange) & (-conf.chunkRange)) * 2;
    }
    if (total_bits < 20) {
        fprintf(stderr, "WARNING: Not enough bits to narrow down structure seed! (%d/%d)\n", total_bits, 20);
        // exit(1);
    }
}

uint64_t gcd(uint64_t a, uint64_t b) {
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    if (a == b)
        return a;
    if (a > b)
        return gcd(a - b, b);
    return gcd(a, b - a);
}

__device__ __host__ void goBack2(uint64_t* rand) {
    *rand = (*rand * 254681119335897ULL + 120305458776662ULL) & MASK48;
}

__device__ __host__ void goBack(uint64_t* rand) {
    *rand = (*rand * 246154705703781ULL + 107048004364969ULL) & MASK48;
}

// stolen from https://github.com/michel-leonard/C-MathSnip/blob/main/mod_inv.c :)
uint64_t modinv(uint64_t ra, uint64_t rb) {
    uint64_t rc, sa = 1, sb = 0, sc, i = 0;
    if (rb > 1) do {
            rc = ra % rb;
            sc = sa - (ra / rb) * sb;
            sa = sb, sb = sc;
            ra = rb, rb = rc;
        } while (++i, rc);
    sa *= (i *= ra == 1) != 0;
    sa += (i & 1) * sb;
    return sa;
}

typedef struct {
    uint64_t *seeds;
    size_t size;
} Seeds;

Seeds solve_constraints_for_lower20(Constraints constraints) {
    int max = 0;
    for (size_t i = 0; i < constraints.length; i++) {
        if ((int)log2(constraints.cons[i].chunk_range & -constraints.cons[i].chunk_range)) {
            max = MAX(max, constraints.cons[i].chunk_range);
        }
    }
    int mp = max & -max; 
    int lp = (int)log2(mp);

    size_t num_lower = 0;
    uint64_t *valid_lower = (uint64_t *)malloc(sizeof(uint64_t) * (int)((100.0/100.0) * (1<<(17 + lp))));
    // huge thank you to Kris for his great writeup on bitlifing: https://github.com/Kludwisz/BitLifting/
    for (uint64_t lower = 0; lower < 1ull<<(17 + lp); lower++) {
        for (size_t i = 0; i < constraints.length; i++) {
            Constraint c = constraints.cons[i];
            if (c.tp != STRUCTURE) {
                continue;
            }
            int p = c.chunk_range & (-c.chunk_range);
            if (p < 2) {
                continue;
            }
            uint64_t lower_seed = ((uint64_t)c.reg_x*341873128712ULL + (uint64_t)c.reg_z*132897987541ULL + lower + (uint64_t)c.salt);
            setSeed(&lower_seed, lower_seed);
            if (nextInt(&lower_seed, c.chunk_range) % p != c.offset_x % p || nextInt(&lower_seed, c.chunk_range) % p != c.offset_z % p) {
                goto next_lower;
            }
        }
        valid_lower[num_lower] = lower;
        num_lower++;
next_lower:
        ;
    }

    return (Seeds){.seeds=valid_lower, .size=num_lower};
}

#define LENGTH(x) (sizeof(x) / sizeof(x[0]))
__device__ __managed__ unsigned long long int checked = 0;
__device__ __managed__ unsigned long long int buffer_size = 0;

__global__ void find_structure_seeds(Constraint best_constraint, Constraints constraints, uint64_t region_seed, uint64_t basis, uint64_t end, uint64_t reduced_mod, uint64_t *buffer) {
    uint64_t input_seed = blockDim.x * blockIdx.x + threadIdx.x;
    if (input_seed % reduced_mod != basis || input_seed > end) {
        return;
    } 
    uint64_t full_region_seed = (input_seed << 20) | region_seed; 
    goBack2(&full_region_seed);
    int x = nextInt(&full_region_seed, best_constraint.chunk_range);
    int z = nextInt(&full_region_seed, best_constraint.chunk_range);
    if (x != best_constraint.offset_x || z != best_constraint.offset_z) { // the second part should never really happen, just a sanity check...
        return;
    }
    goBack2(&full_region_seed);

    full_region_seed ^= 0x5deece66d;
    uint64_t structure_seed = (full_region_seed - best_constraint.reg_x * 341873128712ULL - best_constraint.reg_z * 132897987541ULL - (uint64_t)best_constraint.salt) & MASK48;
    for (size_t i = 0; i < constraints.length; i++) {
        Constraint *c = &constraints.cons[i];

        if (c->tp == PILLAR) {
            uint64_t ss = structure_seed;
            if ((nextLong(&ss) & 65535) != (uint64_t)c->type) {
                return;
            }
        }
        if (c->tp != STRUCTURE) {
            continue;
        }

        uint64_t ss = structure_seed;
        setSeed(&ss, c->reg_x*341873128712ull + c->reg_z*132897987541ull + ss + (uint64_t)c->salt);
        if (nextInt(&ss, c->chunk_range) != c->offset_x || nextInt(&ss, c->chunk_range) != c->offset_z) {
            return;
        }
    }
    // atomicAdd(&checked, 1ull);
    uint64_t r = atomicAdd(&buffer_size, 1ull);
    buffer[r] = structure_seed;
}

#define GPU_ASSERT(code) gpuAssert((code), __FILE__, __LINE__)
inline void gpuAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s (code %d) %s %d\n", cudaGetErrorString(code), code, file, line);
    exit(code);
  }
}

Constraints make_gpu_constraints(Constraints cpu_constraints) {
    Constraint *new_constraints;
    cudaMalloc((void **)&new_constraints, sizeof(Constraint) * cpu_constraints.length);
    cudaMemcpy(new_constraints, cpu_constraints.cons, sizeof(Constraint) * cpu_constraints.length, cudaMemcpyHostToDevice); 
    return (Constraints){.cons=new_constraints, .length=cpu_constraints.length};
}

int main(void) {
    Constraint constraints[] = {
        STRUCTURE_CONSTRAINT(Village, -16 + 16, 736),
        STRUCTURE_CONSTRAINT(Igloo, 272, 608),
    };

    cudaSetDevice(0);

    Constraints cons = (Constraints){.cons=constraints, .length=LENGTH(constraints)};
    parse_constraints(cons);

    Seeds lower20_seeds = solve_constraints_for_lower20(cons);
    printf("Found %ld lower20 seeds\n", lower20_seeds.size);

    int max = 0;
    int best_index = -1;
    for (size_t i = 0; i < cons.length; i++) {
        Constraint *c = &cons.cons[i];
        uint64_t param = c->chunk_range / gcd(8, c->chunk_range);
        if (param > (uint64_t)max) {
            best_index = i;
            max = (int)param;
        }
    }
    Constraint best_constraint = cons.cons[best_index];

    Constraints gpu_cons = make_gpu_constraints(cons);

    FILE *seed_file = fopen("structure_seeds.txt", "w");
    for (size_t i = 0; i < lower20_seeds.size; i++) {
        uint64_t lower20 = lower20_seeds.seeds[i];
        uint64_t region_lower20 = (((uint64_t)best_constraint.reg_x*341873128712ULL + (uint64_t)best_constraint.reg_z*132897987541ULL + lower20 + (uint64_t)best_constraint.salt)) & 0xFFFFF;
        setSeed(&region_lower20, region_lower20);
        int x = nextInt(&region_lower20, best_constraint.chunk_range);
        int z = nextInt(&region_lower20, best_constraint.chunk_range);
        (void)x;
        (void)z;
        region_lower20 = region_lower20 & 0xFFFFF; 
        uint64_t L3 = (region_lower20 >> 17) & 0x7;
        uint64_t C = best_constraint.chunk_range;
        uint64_t Z = best_constraint.offset_z;
        uint64_t d = gcd(8, C);
        uint64_t reduced_mod = C / d;
        uint64_t inv = modinv(8 / d, C / d);
        uint64_t basis = inv * ((Z - L3) / d) % reduced_mod;
        uint64_t end = (1ull<<28ull) - 1ull;

        uint64_t *seed_buffer;
        GPU_ASSERT(cudaMallocManaged((void **)&seed_buffer, sizeof(uint64_t) * 4000000));
        
        find_structure_seeds<<<2097152, 128>>>(best_constraint, gpu_cons, region_lower20, basis, end, reduced_mod, seed_buffer);
        GPU_ASSERT(cudaDeviceSynchronize()); 
        printf("%ld/%ld => %lld buffer_size: %lld\n", i, lower20_seeds.size, checked, buffer_size);

        for (int i = 0; i < buffer_size; i++) {
            fprintf(seed_file, "%lu\n", seed_buffer[i]);
        }

        buffer_size = 0;
        cudaFree(seed_buffer);
    }
    printf("Found: %lld structure seeds!\n", checked);

    fclose(seed_file);
    free(lower20_seeds.seeds);
    return 0;
}