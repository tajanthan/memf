# Memory Efficient Max Flow for Multi-label Submodular MRFs

## Abstract
Multi-label submodular Markov Random Fields (MRFs) have been shown to be solvable using max-flow based on an encoding of the labels proposed by Ishikawa, in which each variable $X_i$ is represented by nodes (where is the number of labels) arranged in a column. However, this method in general requires $2\ell^2$ edges for each pair of neighbouring variables. This makes it inapplicable to realistic problems with many variables and labels, due to excessive memory requirement. In this paper, we introduce a variant of the max-flow algorithm that requires much less storage. Consequently, our algorithm makes it possible to optimally solve multi-label submodular problems involving large numbers of variables and labels on a standard computer.

## Publication
Memory Efficient Max Flow for Multi-label Submodular MRFs.  
Thalaiyasingam Ajanthan, Richard Hartley, and Mathieu Salzmann.  
IEEE Conference on Computer Vision and Pattern Recognition (CVPR), June 2016.  

[pdf][1] [supplementary][2] [poster][3] [code][4]

[1]: docs/memf.pdf "pdf"
[2]: docs/memf_supp.pdf "supplementary"
[3]: docs/memf_poster.pdf "poster"
[4]: https://github.com/tajanthan/memf "code"

