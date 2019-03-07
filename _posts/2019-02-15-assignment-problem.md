---
layout: post
category: code
title: The (Linear Sum) Assignment Problem
implementation: /assets/munkres_assignment_algo.py
---

Consider the problem of placing $n$ individuals into $n$ jobs given a matrix $(a_{ij})$ of the rating of individual $i$ at job $j$. The problem is to find an efficient construction of the optimal assignment. In other words, for any square matrix, pick a sum-maximizing(or minimizing) set of elements that do no share a row or column. In the case where the cost marix is binary, this algorithm also solves the problem of picking a maximal set of nonattacking rooks given such an arrangement on a chess board. Scanning all $n!$ possibilities is a bad option because even placing only 50 individuals would require $3 \cdot 10^{64}$ cases to check. The algorithm described here, sometimes called the _Hungarian Algorithm_ or the _Kuhn-Munkres_ Algorithm, runs in $O(n^3)$ which corresponds to $\approx 10^5$ operations for assigning 50 individuals. Quite an improvement! I first studied this problem as an extremely simple case of a problem in geometric analysis called Monge-Kantorovich Mass Transfer. [A paper](https://pdfs.semanticscholar.org/848c/717ba51e48afef714dfef4bd6ab1cc050dab.pdf) by the great topologist James Munkres outlines a solution, from which I adapted this [Python implementation](/assets/munkres_assignment_algo.py).

These types of problems were the subject of the 1975 Nobel Prize in Economics won by Leonid Kantorovich and Tjalling Koopmans. 
