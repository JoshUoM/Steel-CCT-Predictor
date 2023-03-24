# CCT Predictor for Low Alloy Steels

Python code for predicting the continuous cooling transformation (CCT) behaviour of low alloy steels.

## 1. Inputs

The code requires **3 inputs** in order to run. These are; **comp**, **G**, and **rates**.

### Input 1 - comp

An input for the steel alloy composition.

Alloy composition should be inputted in a Python dictonary format, with curly brackets ({}) used to define the dictionary space. Dictionaries are made up of key:value pairs, where keys and values are separated by a colon (:), and a comma (,) is used to separate key:value pairs. In this case, alloying element symbols are inputted as the keys and their associated concentration in wt.% inputted as the values. 

An example **comp** input would be as follows:

    comp = {'C':0.1,'Si':0.2,'Mn':0.3,'Ni':0.4,'Cr':0.5,'Mo':0.6}
    
for an Fe-0.1C-0.2Si-0.3Mn-0.4Ni-0.5Cr-0.6Mo alloy.

### Input 2 - G

An input for the austenite grain size. 

Austenite grain size should be inputted as the ASTM grain size number. Both a conversion table and equations for calculating this value from SI units can be found in ASTM E112: Standard Test Methods for Determining Average Grain Size.

An example input for **G** would be as follows:

    G = 10
    
for ASTM grain size 10.

### Input 3 - rates

An input for the cooling rates to be tested. 

Cooling rates are to be inputted in a Python list format, with square brackets ([]) used to define the list space. Lists are made up of values separated by a comma (,). In this case, the desired cooling rates should be inputted in degrees per second (°C/s).

An example input for **rates** would be as follows:

    rates = [0.1, 1, 10, 100]
    
for cooling rates 0.1, 1, 10, and 100°C/s.

## 2. Running the Code

The model can be run using a variety of functions. These are 'Steel_CCT_Calculator', for predicting the CCT curves, 'CCT_Plotter', for plotting the CCT curves from the first function, and 'CCT_Fractions', for calculting the final constituent fractions at each modelled cooling rate.

### Function 1: CCT Calculator

This function can be used to get the raw CCT data for the test, and is as follows:

    Steel_CCT_Calculator(comp,G,rates)

Example code:

    comp = {'C':0.1,'Si':0.2,'Mn':0.3,'Ni':0.4,'Cr':0.5,'Mo':0.6}
    G = 10
    rates = [0.1, 1, 10, 100]
    
    Steel_CCT_Calculator(comp,G,rates)
    
### Function 2: CCT Plotter

This function can be used to plot the raw data as a CCT curve (as shown in Fig.1). The function is as follows:

    CCT_Plotter(Ts,comp,rates)
    
where **Ts** is the direct output taken from the function 'Steel_CCT_Calculator'.

Example input:
    
    comp = {'C':0.1,'Si':0.2,'Mn':0.3,'Ni':0.4,'Cr':0.5,'Mo':0.6}
    G = 10
    rates = [0.1, 1, 10, 100]
    
    Ts = Steel_CCT_Calculator(comp,G,rates)
    
    CCT_Plotter(Ts,comp,rates)
    
Example output:

<figure>
  <img
  src="example figures/CCT_Plotter_EXAMPLE.png"
  alt="."
  width="50%" 
  height="50%">
  <figcaption>Fig. 1 An example output when using the 'CCT_Plotter' function. Constituent nomenclature is as follows: f - ferrite, p - pearlite, b - bainite, m - martensite.</figcaption>
</figure>

&nbsp;
   
### Function 3: CCT Fractions

This function can be used to output the final constituent fractions, X, from the modelled CCT curves (as shown in Fig. 2) and is used as follows:

    CCT_Fractions(Ts,rates)
    
where **Ts** is the direct output from function 'Steel_CCT_Calculator'.

Example input:

    comp = {'C':0.1,'Si':0.2,'Mn':0.3,'Ni':0.4,'Cr':0.5,'Mo':0.6}
    G = 10
    rates = [0.1, 1, 10, 100]
    
    Ts = Steel_CCT_Calculator(comp,G,rates)
    
    CCT_Fractions(Ts,rates)
    
Example output:

<figure>
  <img
  src="example figures/CCT_Fractions_EXAMPLE.png"
  alt="."
  width="50%" 
  height="50%">
  <figcaption>Fig. 2 An example output when using the 'CCT_Fractions' function. Output are the final fractions, X, for each constituent. Constituent nomenclature is as follows: f - ferrite, p - pearlite, b - bainite, m - martensite, a - austenite.</figcaption>
</figure>

&nbsp;
