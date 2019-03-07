---
layout: post
category: code
title: Computing The Convex Hull in $\mathbb{R}^2$ and $\mathbb{R}^3$
implementation: /assets/munkres_assignment_algo.py
---

There are many classical algorithms to compute the convex hull of a finite point set in the plane, but I chose to implement three. Jarvis Matching involves iteratively adding points to the hull as to mimimize a certain turning angle. It runs in $O(nh)$ where $h$ is the number of points in the hull. The Graham scan iterates through groupings of three points naturally induced by modified polar coordinates, allowing it to reduce to two simple cases. Ronald Graham's algorithm runs in $O(n \log(n))$. 