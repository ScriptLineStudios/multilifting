# multilifing
implementation of multilifting for Minecraft structure cracking 

# what is multilifting?
multilifting is a technique which combines bit lifting and state prediction to crack strcuture seeds from a set of constraints, using a regular bitlifting approach will generally result in set of the lower 20 bits of structure seeds, from there it may like a 28 bit bruteforce is needed for each lower20 in order for a full structure seed to be found, however we can do better, depending on the structures constraining the search, the 2^28 search space can be reduced as follows:
| Structure     | Search space |
| ------------  | ------------ |
| Ruined Portal | ~2<sup>23.35</sup>     |       
| Desert Temple | ~2^26.41     |           
| Jungle Temple | ~2^26.41     |         
| Shipwreck     | ~2^25.67     |         
| Village       | ~2^24.29     |         
