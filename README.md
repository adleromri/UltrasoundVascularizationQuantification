# UltrasoundVascularizationQuantification
GUI tool for analyzing and quantifying the vascularization within a given ultrasonic video clip taken with both normal B-Mode and harmonic imaging at the same time.

PanGUI_1_4(new).mlapp - The main GUI file.
PanGUI_TGC.m - Time Gain Compensation function (Gamma correction).
PanGUI_RingSegImg.m - Segments the image into its given rings.
PanGUI_CC.m - Cross-Correlation function (frames selection).
PanGUI_Calc.m - Calculates the rings' values and temporal properties.
PanGUI_RingsPlot.m - Creates .jpg files from the rings images.
PanGUI_Plot_2.m - Creates .jpg files from the rings' analysis images.
PanGUI_Reg.m - Registration function.
PanGUI_Dec.m - Deconvolution main function. (Calls: "ker_check_2.m" and "decon5_p3.m")
ker_check_2.m - Kernel calculation function. (Calls: "AlgFunc5.m")
AlgFunc5.m - Deconvolution image calculation function.
decon5_p3.m - Deconvolution video calculation function. (Calls: "AlgFunc5_2.m")
AlgFunc5_2.m - Deconvolution image calculation function.
