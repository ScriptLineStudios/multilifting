# multilifing
implementation of multilifting for Minecraft structure cracking 

# what is multilifting?
multilifting is a technique which combines bit lifting and state prediction to crack strcuture seeds from a set of constraints, using a regular bitlifting approach will generally result in set of the lower 20 bits of structure seeds, from there it may seem like a 28 bit bruteforce is needed for each lower20 in order for a full structure seed to be found, however we can do better, depending on the structures constraining the search, the 2<sup>28</sup> search space can be reduced as follows:
| Structure     | Search space | | Reduction
| ------------  | ------------ | |
| Ruined Portal | ~2<sup>23.35</sup>     | |
| Village       | ~2<sup>24.29</sup>     | | 
| Trail Ruins | ~2<sup>24.29</sup>      |  |
| Shipwreck     | ~2<sup>25.67</sup>     |  |       
| Swamp Hut | ~2<sup>26.41</sup>      |   |
| Desert Temple | ~2<sup>26.41</sup>      | |           
| Jungle Temple | ~2<sup>26.41</sup>     |  | 
| Igloo | ~2<sup>26.41</sup>      |   |
