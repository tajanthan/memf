# Memory Efficient Max Flow

## Related Publication

This code implements the MEMF algorithm	described in the following paper  

[Memory Efficient Max Flow for Multi-label Submodular MRFs](https://tajanthan.github.io/memf).  
Thalaiyasingam Ajanthan, Richard Hartley, and Mathieu Salzmann.  
IEEE Conference on Computer Vision and Pattern Recognition (CVPR), June 2016.  
IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), March 2018.

## Code Assumptions

1. The MRF energy has the following form  
   $$E(x) = \sum \theta_{i}(x_i) + \sum \theta_{ij} (x_i, x_j)\ ,$$  
   where $\theta_{ij} (x_i, x_j) = \gamma_{ij} \theta(|x_i - x_j|)$.
2. The code currently supports MRF with 4-connected grid structure only and nodes are labelled from 0 --> width * height - 1 in    raw major ordering.  


## Contact

This code is for research purposes only, if you want to use it for commercial purpose please contact us.  
**Email:** ajanthan {at} robots {dot} ox {dot} ac {dot} uk

## Example Usage

To assist the user, example.cpp, Makefile and sample data files are provided. The following command runs the MEMF algorithm on 10x10 image with 5 labels, with quadratic pairwise potential.

``` 
memf.exe 10 10 5 <sample>/toy_unary_10_10_5.txt <sample>/toy_binary_4_10_10.txt <sample>/toy_binaryPot_10_10_5_l2.txt
```
[Code Ocean Link](https://codeocean.com/2018/03/22/memf-colon-memory-efficient-max-flow-for-multi-label-submodular-mrfs/code)
