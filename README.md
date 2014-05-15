How to use it: 
===============

* Content :
  - JPsiEvents.C will create from the analyzer output a miniTree with 1 entry by  GEN muon and containing matching info with STA/seed + in case of matching, basic info on the seed
  - createEfficiencyPlotsJPsi.C will create the efficiency plots from the miniTree
  - example_of_analysisOutput.root is a example of analyzer output 

* How to use it ?
  - `root -b -q -l 'JPsiEvents.C("example_of_analysisOutput.root","exampleOfMinitree")'` -> it will create a file called minitree_exampleOfMinitree.root 
  - `root -l -b -q 'createEfficiencyPlotsJPsi.C("exampleOfMinitree")'` -> will create a file histo_exampleOfMinitree.root containing all the efficiency plots and some histograms
