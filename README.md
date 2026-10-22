# Neurophysiological-Correlates-of-Cortical-Hierarchy-across-the-Lifespan

This repository contains code and data for processing, statistical analyses and visualization in support of the work "Neurophysiological correlates of cortical hierarchy across the lifespan" 
Additionally, Code for preprocessing of MEG and MRI data, cortical thickness estimation as well as most of the time-series feature computation can be found here: https://github.com/chstier/Aging_timefeats.

## Data availability
Raw data can be obtained via: https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/ upon request and under specified conditions.

## Toolboxes
The complete analysis was performed using MATLAB R2018b, R2022b, Python3.11, RStudio 2023.03.0:

Most time-series features were computed using the hctsa toolbox (https://hctsa-users.gitbook.io/hctsa-manual). To obtain a parcellated map of cortical hierarchy and perform spin testing, we used the neuomaps toolbox (https://netneurolab.github.io/neuromaps/). Visualization of features values along the human cortex were performed using the brainstat toolbox  in Matlab R2022b 94 and the scientific color maps 8.0.1.

## References
- B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
- B.D. Fulcher, M.A. Little, N.S. Jones. Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 20130048 (2013).
- Markello, R. D. et al. neuromaps: structural and functional interpretation of brain maps. Nat. Methods 19, 1472–1479 (2022).
- Crameri, F., Shephard, G. E. & Heron, P. J. The misuse of colour in science communication. Nat. Commun. 11, 5444 (2020).
