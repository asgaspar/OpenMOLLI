#   <div align="center">  Open-MOLLI </div>
##  <div align="center"> Open-source Inversion recovery myocardial T1 mapping sequence

<div align="justify"> Quantitative MR methods require reproducibility studies to evaluate their accuracy and precision which can be difficult when MR sequences can vary between centers and scanners. In addition, faster prototyping for concept testing and improvement is challenging and time consuming.</div>

<br/>

This package offers an open-source Prototype of Myocardial T1 mapping (Open-MOLLI) using Pulseq [[1,2]](#references) which includes an inversion recovery T1 mapping sequence with a triggering scheme. Open-MOLLI accessibility allows faster implementation of new ideas, while its applicability to different vendors though Pulseq versatility makes it an easier route for reproducibility studies. 


<p align="center">
<img src="OpenMOLLI_arial.png"/>
</p>

<br/>

If you use the sequence Open-MOLLI in your work, cite as:

```
Gaspar AS, Silva NA, Ferreira AM, Nunes RG. Repeatability of Open-MOLLI: An open-source inversion recovery myocardial T1 mapping sequence for fast prototyping. Magn Reson Med. 2024; 1-10. doi: 10.1002/mrm.30080.
```

## Packages
Open-MOLLI can be build with Matlab or Python: 
*  **Matlab**: 
	* **Open-MOLLI**  can be build from code in Matlab_Open-MOLLI folder. Files in mr+ folder should be added Pulseq mr+ folder. 
* **Python**:  
	* **PyOpen-MOLLI** can be build from code in Python_pyOpen-MOLLI folder. This requires pypulseq [[2]](#references) . 
	* A tutorial notebook for pyOpen-MOLLI is available at [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/asgaspar/OpenMOLLI/blob/main/PyOpenMOLLI_Python/pyOpenMOLLI.ipynb)


## Requirements
In order to create a `Open-MOLLI.seq` file you will need: 
*  **Matlab**:  **Pulseq** package available at: https://github.com/pulseq/pulseq
*   **Python**:  
	* **pyPulseq** [[2]](#references)  package available at: https://github.com/imr-framework/pypulseq 
	* You can also install with `pip install pypulseq`

## Utilities Reconstruction
In order to use the reconstruction utilities from repository https://github.com/fluese/reconstructionPipeline/
*  Read *.dat:  **mapVBVD** package available at: https://github.com/fluese/reconstructionPipeline/tree/master/externalTools/mapVBVD
*  GRAPPA:   **GRAPPA** available at: https://github.com/fluese/reconstructionPipeline/tree/master/externalTools/GRAPPA
*  Filter:   **Tukey** available at: https://github.com/fluese/reconstructionPipeline/blob/master/processingPipeline/window3_tukeytaper.m


## References
1. Layton KJ, Kroboth S, Jia F, Littin S, Yu H, Leupold J, Nielsen JF, Stöcker T and Zaitsev M. Pulseq: A rapid and hardware‐independent pulse sequence prototyping framework. Magn Reson Med. 2017;77:1544-1552. https://doi.org/10.1002/mrm.26235
2. Keerthi R, Geethanath S, and Vaughan J. PyPulseq: A Python Package for MRI Pulse Sequence Design. Journal of Open Source Software. 2019;4(42): 1725. https://doi.org/10.21105/joss.01725
