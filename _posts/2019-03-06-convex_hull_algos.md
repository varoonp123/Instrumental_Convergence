---
layout: post
category: code
title: Computing The Convex Hull in $\mathbb{R}^2$
implementation: /assets/munkres_assignment_algo.py
---

There are many classical algorithms to compute the convex hull of a finite point set in the plane, but I chose to implement three. Jarvis Matching involves iteratively adding points to the hull as to mimimize a turning angle. It runs in $O(nh)$ where $h$ is the number of points in the hull. The [Graham scan](http://www.math.ucsd.edu/~ronspubs/72_10_convex_hull.pdf) iterates through groupings of three points naturally induced by modified polar coordinates, allowing it to reduce to two simple cases. Ronald Graham's algorithm runs in $O(n \log(n))$. I created these [Python implementations](/assets/jarvismatch.py) from the original papers. 

<div class="convexhull_pics">
    <figure>
        <img src="/assets/images/convexhull_7pts_scaled.jpg" alt="The outline of the convex hull of seven planar points."/>
        <figcaption>The outline of the convex hull of seven planar points.</figcaption>
    </figure>
    <figure>
        <img src= "/assets/images/convexhull_250kpts_scaled.jpg" alt="Convex hull of 250,000 points."/> 
        <figcaption>Convex hull of 250,000 points.</figcaption>
    </figure>
</div>

