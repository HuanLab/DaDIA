**Data-Dependent-Assisted Data-Independent Acquisition (DaDIA.R) User
Manual**

(Version 4.0; Sep 21^th^, 2020)

Jian Guo, Sam Shen, Tao Huan

-   DaDIA.R is a metabolic feature extraction and annotation workflow
    for DaDIA metabolomics workflow.

-   The program is written in the language 'R' and is publicly available
    at GitHub ([https://github.com/HuanLab/DaDIA.git]{.ul})

-   Please see below for instructions on using the DaDIA.R code:

1)  Download the R-scrip "DaDIA_SWATHprocessing.R" from Github

2)  Install libraries "XCMS", "MSnbase", "dplyr", "doParallel",
    "foreach", "metaMS", and "CAMERA" if previously not installed (R
    Version 4.0 or above and XCMS Version 3.11.4 is required; all other
    packages should be updated to the newest available version)

3)  In line 17 -- 63, set required parameters to desired value according
    to onscreen prompts (See below).

![A screenshot of a social media post Description automatically
generated](media/image1.png){width="6.5in" height="6.6875in"}

**Table 1.** DaDIA parameters needed to be changed according to the
user\'s preference.

  **Row \#**   **Parameter Name**            **Parameter Function**
  ------------ ----------------------------- -----------------------------------------------------------------------------------------------------------------------------------
  18           *DDA.directory*               Set directory containing all DDA .mzxml files
  19           *DIA.directory*               Set directory containing all DIA(SWATH) .mzxml files, SWATH isolation window labeling .txt file, and annotation library .msp file
  20           *cwpDDA*                      Set XCMS parameters for DDA feature extraction
  27           *cwpDIA*                      Set XCMS parameters for DIA feature extraction
  34           *mass.tol*                    Set m/z tolerance (± ppm) for feature dereplication and MS2 matching
  35           *mass.const.tol*              Set m/z tolerance (± constant value) for feature rescue
  36           *rt.tol*                      Set retention time tolerance (± sec) for identifying same features
  37           *num.samples*                 Set number of DIA(SWATH) samples to run
  38           *plot.DaDIA*                  Set whether to plot extracted DaDIA MS1 features
  39           *plot.DaDIA.mztol*            Set DaDIA feature plotting mz window width
  40           *plot.DaDIA.rttol*            Set DaDIA feature plotting rt window width
  42           *bw*                          Set XCMS feature alignment bandwidth
  43           *minfrac*                     Set XCMS feature alignment minimum sample fraction
  44           *mzwid*                       Set XCMS feature alignment m/z slice width
  45           *max*                         Set XCMS feature alignment maximum \# of groups / slice
  46           *quantitative.method*         Set whether to use max or integrated intensity for calculations; pick between peak height or peak area
  51           *feature.annotation*          Set whether to perform MS2 extraction and DaDIA feature annotation
  52           *db.name*                     Set the name of the library used for dot product annotation
  53           *RDS*                         Set whether the annotation library is in RDS of MSP format
  56           *ms1.tol*                     Set MS1 tolerance in dot product calculation
  57           *ms2.tol*                     Set MS2 tolerance in dot product calculation
  58           *dot.product.threshold*       Set annotation dot product score threshold
  59           *match.number.threshold*      Set annotation match number score threshold
  60           *adduct_isotope.annotation*   Set whether to perform CAMERA adduct and isotope annotation
  61           *export.mgf*                  Set whether to export MS2 spectrum as individual .mgf files
  62           *combine.mgf*                 Set whether to concatenate all exported .mgf files into a single .mfg file
  63           *annotation.plot*             Set whether to plot all annotated MS2 spectrum against library

4)  In line 7, set the directory in the user's computer that contains
    all DDA samples (Note: there should be only .mzxml files in this
    folder). In line 8, set the directory in the user's computer that
    contains all SWATH samples, SWATH pocket definitions in .txt format,
    and the annotation library in .msp format (Note: there should only
    be .mzxml, .txt, and .msp files in this folder). See below for a
    sample format of the SWATH pocket definition file; note that the
    column headers should be kept the same as the example shown:

![A screenshot of a cell phone Description automatically
generated](media/image2.png){width="2.2336450131233594in"
height="2.9232261592300963in"}

5)  Note: if users wish to use MS-Dial library for annotation, they must
    first convert the raw MS-Dial library to a readable .msp format
    using the included R script *"convertMSP.R"*

6)  Click on "Source" in the R-Studio interface to begin the DaDIA
    workflow

7)  After running the scrip for single sample runs, a csv file
    "DaDIAtable.csv" containing all MS1 level features extracted and a
    csv file "annotated_output.csv" containing feature annotation
    results will be generated. After running the scrip for multi-sample
    runs, multiple csv files "n_DaDIAtable.csv" (n is number of
    DIA(SWATH) samples) containing all MS1 level features extracted for
    each DIA(SWATH) sample, a csv file "alignedDaDIAtable.csv"
    containing all MS1-features in different samples aligned together,
    and a csv file "annotated_output.csv" that adds annotation to the
    "alignedDaDIAtable" will be generated. If the user wishes to plot
    DaDIA features or output .mgf files, additional folders and files
    will be generated. All results generated by the program will be in
    the *DIA.directory* set by the user.
