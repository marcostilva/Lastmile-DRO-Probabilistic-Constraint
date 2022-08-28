# Lastmile-DRO-Probabilistic-Constraint
Stochastic crowd shipping last-mile delivery with correlated marginals and probabilistic constraints

In this work, we study last-mile delivery with the option of crowd shipping.
A company uses occasional drivers to complement its fleet in the activity of
delivering products to its customers. We model it as a variant of the stochastic
capacitated vehicle routing problem. Our approach is data-driven, where not
only customer orders but also the availability of occasional drivers are uncertain.
It is assumed that marginal distributions of the uncertainty vector are known,
but the joint distribution is difficult to estimate. We optimize considering a
worst-case joint distribution and model with a strategic planning perspective,
where we calculate an optimal a priori solution before the uncertainty is revealed.
A limit on the infeasibility of the routes due to the capacity is imposed using
probabilistic constraints.
We propose an extended formulation for the problem using column-dependent
rows and implement a branch-price-and-cut algorithm to solve it. We also develop 
a heuristic approximation to cope with larger instances of the problem.
