# SPTAnalyser

![](tmp/SPTAnalyser_software_workflow.png)

Batch processing of single particle tracking data</br>
Compatible with PALMTracer<sup>[1]</sup>, rapidSTORM<sup>[2]</sup>, ThunderSTORM<sup>[3]</sup> and swift<sup>[4]</sup></br>
MSD-based extraction of diffusion coefficients and diffusion types<sup>[5]</sup></br>
Transition counting between different diffusion modes</br>
Hidden Markov modeling with ermine<sup>[6]</sup>

## To install the package run:

```
conda create --name SPTAnalyser
conda activate SPTAnalyser
conda install pip pywin32
pip install jupyterlab
python -m pip install pyErmine
pip install SPTAnalyser-XXX-py3-none-any.whl
pip install jupyter_contrib_nbextensions
jupyter notebook
```

Give the analysis a try with the test files in the dataset folder, including localization files from ThunderSTORM, tracked files from swift, and SPTAnalyser trackAnalysis output.

## Contributors
Johanna Rahm, Sebastian Malkusch, Marie-Lena Harwardt, Marina Dietz, Claudia Catapano, Alexander Niedrig

## Literature

[1] https://docs.google.com/forms/d/e/1FAIpQLSdxkJXmoc5uRqoL8kzCSV2Rv90kEfRMIoaPuk--Bpm0Lttb0g/viewform </br>
[2] S. Wolter, A. Löschberger, T. Holm, S. Aufmkolk, M.-C. Dabauvalle, S. van de Linde, M. Sauer, 2012, Nature Methods, 9, 1040-1041, DOI: 10.1038/nmeth.2224 </br>
[3] M. Ovesny, P. Krizek, J. Borkovec, Z. Svindrych, G. M. Hagen, 2014, Bioinformatics, 30, 2389-2390, DOI: 10.1093/bioinformatics/btu202 </br>
[4] M. Endesfelder, C. Schießl, B. Turkowyd, T. Lechner, U. Endesfelder, Manuscript in Prep.; http://bit.ly/swifttracking </br>
[5] M.-L. I. E. Harwardt, P. Young, W. M. Bleymüller, T. Meyer, C. Karathanasis, H. H. Niemann, M. Heilemann, M. S. Dietz, 2017, FEBS Open Bio, 7, 1422-1440, DOI: 10.1002/2211-5463.12285 </br>
[6] https://github.com/SMLMS/pyErmine

## Citation
Please cite our paper introducing the workflow, if you use SPTAnalyser for your research. </br>
J. V. Rahm, S. Malkusch, U. Endesfelder, M. S. Dietz, M. Heilemann, Front. Comput. Sci., 12 November 2021, https://doi.org/10.3389/fcomp.2021.757653

## Web
For more information about the authors visit:  </br>
http://share.smb.uni-frankfurt.de  </br>

More single-particle tracking data can be found at: </br>
https://www.ebi.ac.uk/biostudies/studies/S-BSST712 </br>
