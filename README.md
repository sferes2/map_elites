# map_elites

MAP-Elites module for Sferes2.

## Authors:
- Jean-Baptiste Mouret jean-baptiste.mouret@inria.fr
- Danesh Tarapore

## Reference:

Mouret, J. B., & Clune, J. (2015). Illuminating search spaces by mapping elites. arXiv preprint arXiv:1504.04909.

http://arxiv.org/abs/1504.04909

Many fields use search algorithms, which automatically explore a search space to find high-performing solutions: chemists search through the space of molecules to discover new drugs; engineers search for stronger, cheaper, safer designs, scientists search for models that best explain data, etc. The goal of search algorithms has traditionally been to return the single highest-performing solution in a search space. 

Here we describe a new, fundamentally different type of algorithm that is more useful because it provides a holistic view of how high-performing solutions are distributed throughout a search space. It creates a map of high-performing solutions at each point in a space defined by dimensions of variation that a user gets to choose. This Multi-dimensional Archive of Phenotypic Elites (MAP-Elites) algorithm illuminates search spaces, allowing researchers to understand how interesting attributes of solutions combine to affect performance, either positively or, equally of interest, negatively. 

For example, a drug company may wish to understand how performance changes as the size of molecules and their cost-to-produce vary. MAP-Elites produces a large diversity of high-performing, yet qualitatively different solutions, which can be more helpful than a single, high-performing solution. Interestingly, because MAP-Elites explores more of the search space, it also tends to find a better overall solution than state-of-the-art search algorithms. We demonstrate the benefits of this new algorithm in three different problem domains ranging from producing modular neural networks to designing simulated and real soft robots. 

Because MAP- Elites (1) illuminates the relationship between performance and dimensions of interest in solutions, (2) returns a set of high-performing, yet diverse solutions, and (3) improves finding a single, best solution, it will advance science and engineering. 