---
layout: post
category: code
title: The (Linear Sum) Assignment Problem
implementation: /assets/munkres_assignment_algo.py
---

Consider the problem of placing $n$ individuals into $n$ jobs given a matrix $(a_{ij})$ of the rating of individual $i$ at job $j$. The problem is to find an efficient construction of the optimal assignment. In other words, for any square matrix, pick a sum-maximizing(or minimizing) set of elements that do no share a row or column. In the case where the cost marix is binary, this algorithm also solves the problem of picking a maximal set of nonattacking rooks given such an arrangement on a chess board. I created a [Python implementation](/assets/munkres_assignment_algo.py) of the _Hungarian Algorithm_ or the _Kuhn-Munkres Algorthm_ from a [paper](https://pdfs.semanticscholar.org/848c/717ba51e48afef714dfef4bd6ab1cc050dab.pdf) by James Munkres. These types of problems were the subject of the 1975 Nobel Prize in Economics won by Leonid Kantorovich and Tjalling Koopmans. 


Examples: [1](/assets/assignment_small.txt), [2](/assets/assignment_mid.txt), [3](/assets/assignment_large.txt)

Solutions: [1](/assets/assignment_small_solution.txt), [2](/assets/assignment_mid_solution.txt), [3](/assets/assignment_large_solution.txt)


