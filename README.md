# WallE
Bounding Box Optimizer for large problems where the optimal solution lies on the boundary.

This algorith is a very simple Steepest Descent with projection on the Bounding Box, with some modifications 
like the use of Conjugate Gradients whenever possible (the set of blocked variables do not change in subsequent
iterations) and adaptable moving limits. The line-search is a very crude simple search and accepts the solution 
whenever the objective function decreases (we are not testing for Wolf's condition). 

I will post the theory behind the algorithm soon.
