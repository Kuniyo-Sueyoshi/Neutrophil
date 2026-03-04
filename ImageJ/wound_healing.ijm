run("8-bit");
setAutoThreshold("Default dark no-reset");
//run("Threshold...");
resetThreshold;
run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 autoscale saturate");
setAutoThreshold("Default dark no-reset");
//run("Threshold...");
setThreshold(75, 255, "raw");
setThreshold(75, 255, "raw");
setThreshold(75, 150, "raw");
setThreshold(75, 150, "raw");
//setThreshold(75, 150);r
setOption("BlackBackground", true);
run("Convert to Mask");
run("Minimum...", "radius=6");
