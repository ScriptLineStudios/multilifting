# multilifing
implementation of multilifting for Minecraft structure cracking 

# what is multilifting?
multilifting is a technique which combines bit lifting and state prediction to crack structure seeds from a set of constraints, using a regular bitlifting approach will generally result in set of the lower 20 bits of structure seeds, from there it may seem like a 28 bit bruteforce is needed for each lower20 in order for a full structure seed to be found, however we can do better, depending on the structures constraining the search, the 2<sup>28</sup> search space can be reduced as follows:
| Structure     | Search space       | Reduction |
| ------------  | ------------------ | --------- |      
| Ruined Portal | ~2<sup>23.35</sup> | 25        |
| Village       | ~2<sup>24.29</sup> | 13        |
| Trail Ruins   | ~2<sup>24.29</sup> | 13        |
| Shipwreck     | ~2<sup>25.67</sup> | 5         |   
| Swamp Hut     | ~2<sup>26.41</sup> | 3         |
| Desert Temple | ~2<sup>26.41</sup> | 3         |
| Jungle Temple | ~2<sup>26.41</sup> | 3         |
| Igloo         | ~2<sup>26.41</sup> | 3         |
