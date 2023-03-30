# Athena Analysis
A Python code base for analyzing Athena++ data efficiently 

**Install Instructions**
 1. Clone the repo to wherever you'd like it.
 2. Configure file_config.py, you will find in the roots directory (which is inside the branches directory). There are instructions in the file on what to fill out. Note Athena Analysis expects you to work with multiple datasets where each data set is identified by a string known as its dname, you will pass functions the dname when you call them which will tell Athena_Analysis where to find the data and where so save the results etc.
 3. Configure params.py (you will find it in roots with file_config.py) to match your simulation's parameters. This code base was originally written with a compact binary simulation in mind so if your simulation differs significantly you may not need all the parameters listed, additionally you may want to add your own parameters for aspects of your simulation not covered by the parameters already there. While the Athena_Analysis code base won't call any of these additional parameters on its own you can still call and use them in your own code, so for organization purposes we recommend storing all physical parameters of your simulation here (note whenever you call the main athena_analysis root, or any of the branches a Param object called sim will be created so any parameter in the params file can be accessed as an attribute of sim ie sim.my_parameter)
 4. You're good to go you may now code in leaf.py or other files in the same directory. leaf.py is meant to show you where in the file hierarchy your code should go, you do not actually have to code in the leaf.py file.