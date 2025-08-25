# CNM
Depression core network model (CNM)

The CNM seed regions and the masks for dorsolateral prefrontal cortex (DLPFC) mask are available publicly in NeuroVault (https://neurovault.org/collections/21676/).

If you use the data from this collection please include the following persistent identifier in the text of your manuscript: https://identifiers.org/neurovault.collection:21676. This will help to track the use of this data in the literature.

**The code in this repository was used for producing the results in:**

Komulainen E, Salminen A-S, Aydogan DB, Pamilo S, Raij TT. Cingulo-opercular connectivity enhances the repeatability of transcranial magnetic stimulation target map. (*under review*)


## Instructions for use

Code is shared under the `code` folder of the repository. The python codes can be executed in a virtual environment of Python with necessary packages/modules.

**connectivities_and_correlations.py:**

- Computes mean connectivity maps (target maps) and within-subject spatial correlations between sessions wihtin the DLPFC maks

- Preprocessing and computation of functional connectvity maps of each seed region was computed with DPABI pipeline (rfmri.org)

- Resulting functional connectivity maps have to be arranged: subjectX/Results/FC_FunImg…./FCmaps, subjectX/S2_Results/FC_FunImg…/FCmaps etc.


*e.g. to compute CNM connectivity map:*

```
    from connectivities_and_correlations import compute_mean_connectivity_maps

    # mask.nii=DLPFC-maks
    # ROI1=left dACC
    # ROI2=right dACC
    # ROI3=left AI
    # ROI4=right AI
    # ROI5=left amg
    # ROI6=right amg
    # ROI7=left sgACC
    # ROI8=right sgACC

    . . .

    compute_mean_connectivity_maps(dname='/Path/to/data/', subj_name='subject', include=['ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI7','ROI8'], multiply=['ROI5','ROI6','ROI7','ROI8'], mask_file='/path/to/mask/mask.nii', nvox=100)   

```

**top_voxels_masked.py:**

- Finds coordinates with highest connectivity (i.e., correlation) values (peak coordinates) in each hemisphere, within the DLPFC mask

- Functional connectivity maps need to be at: /path/to/data/subjectx/FCmaps

- DLPFC-mask need to be at: /path/to/data/subjectX/Mask

To execute:
```
    top_voxels_masked.py /path/to/data/
```


**runClusterAnalyses.m:**

- Clusters sgACC and computes connectivity maps from each cluster weighted by cluster size

To execute (with four imaging sessions):

```
    subject='X';

    runClusterAnalysis(subject,'');

    runClusterAnalysis(subject,'S2');

    runClusterAnalysis(subject,'S3');

    runClusterAnalysis(subject,'S4');
```
