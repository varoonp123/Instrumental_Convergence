---
layout: post
category: code
title: The Assignment Problem
implementation: /assets/munkres_assignment_algo.py
---

Consider placing $n$ individuals in $n$ jobs given ratings for each. This gives a matrix $(a_{ij})$ of the rating of individual $i$ at job $j$. Find an efficient method of assigning one person to each job to maximize the total rating. Scanning all $n!$ possibilities is a bad option because even placing only 50 individuals would require $3 \cdot 10^{64}$ cases to check. The algorithm described here runs in $O(n^3)$ which corresponds to $\approx 10^5$ for 50 individuals. Quite an improvement! I first studied this problem as an extremely simple case of a problem in geometric analysis called Monge-Kantorovich Mass Transfer. [A paper](https://pdfs.semanticscholar.org/848c/717ba51e48afef714dfef4bd6ab1cc050dab.pdf) by the great topologist James Munkres outlines a solution, from which I adapted this [Python implementation](/assets/munkres_assignment_algo.py).