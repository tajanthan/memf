# memf
Memory Efficient Max Flow

-------------------------------------------------------------------------------------------
This code implements the MEMF algorithm	described in the following paper 

	"Memory Efficient Max Flow for Multi-label Submodular MRFs", 
	Thalaiyasingam Ajanthan, Richard Hartley and Mathieu Salzmann,
	IEEE Conference on Computer Vision and Pattern Recognition,
	June 2016.

	Code Assumptions:
		1. The MRF energy has the following form
			E(x) = \sum \theta_{i}(x_i) + \sum \theta_{ij} (x_i, x_j),
			where \theta_{ij} (x_i, x_j) = \gamma_{ij} g(x_i, x_j)
		2. The code currently supports MRF with 4-connected grid structure only
			nodes labelled from 0 --> width * height - 1, in a grid structure.
			E.g. width = 3, height = 2
			0 -- 1 -- 2
			|	 |	  |
			3 -- 4 -- 5

	If you use this code, please consider citing the aforementioned paper 
	in any resulting publication.

	This code is for research purposes only, if you want to use it for commercial purpose
	please contact us.
	
	Contact: thalaiyasingam.ajanthan@data61.csiro.au
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
To assist the user, example.cpp, Makefile and sample data files are provided.

	Example usage 
	memf.exe 10 10 5 <sample>/toy_unary_10_10_5.txt <sample>/toy_binary_4_10_10.txt ...
			<sample>/toy_binaryPot_10_10_5_l2.txt
		runs the MEMF algorithm on 10x10 image with 5 labels, with 
		quadratic pairwise potential.
-------------------------------------------------------------------------------------------