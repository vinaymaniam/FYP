There are two versions of U-FRESH code: a fast version and a slow version. 
The algorithm in the ICASSP paper is corresponding to slow version. 
But actually, the fast version can achieve similar results, but with faster speed. 
The difference is that the fast version removes the relibale mapping selection algorithm, but add a second stage upscaling. 
For image compression, if you consider about the speed, I suggest you use the fast version, otherwise the slow version is OK.