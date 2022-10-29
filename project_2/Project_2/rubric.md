# Grading rubric

There are a total of 115 points.

## Code (35 pts)

### Code Style (15 pts)

Each item is graded on a 3-point scale:  
**3** Standard is consistently followed  
**2** Standard is mostly followed, but there are gaps  
**1** Standard is not followed; no consistency  
**0** Standard is not applicable for this project

1. Source code lines are < 80 characters in length; long-lines are broken
2. Blank lines are used to enhance readability (but are not overdone)
3. Code is structured clearly; deeply nested loops and lengthy functions are avoided, where possible
4. Variable and function names are clear, and a consistent style, e.g., "lower\_with\_under" is followed. Constants have a distinguishing style.
5. All functions have a descriptive docstring: arguments and return values for each function are clearly specified.

### Modularity (10 pts)

**10** The project is modular, with common tasks are put into reusable functions and lack or repetition  
**6** Duplication of existing code: for example, constants are (re-)defined in multiple places  
**2** Stringy code with little organization evident

### Efficiency (10 pts)

Code makes use of routines that are tuned for numerical calculation: for example, `numpy` arrays, as opposed to inefficiently looping over python `lists`; iterators are preferred to `while` loops; numbers are handled with scientific notation, rather than as m*10**e.

## Calculation (60 pts)

### Methods (30 pts)

**5 pts** for each item

1. `test_eos.py` reports success
2. Derivatives are correct
3. Central values are correct
4. Lengthscales are correct
5. Integration loop correctly adapts stepsize
6. Pressure guess is correct


### Convergence (15 pts)

**5 pts** for each item

1. &delta; is appropriate&mdash;converged, but not unnecessarily small
2. &eta; is appropriate
3. &xi; is appropriate

**BONUS:** +5 points for calibrated guess for central pressure


### Construction of mass-radius relation (15pts)

**5 pts** for each item

1. Rootfind method is used correctly
2. Central pressure is 0.77 GM<sup>2</sup>/R<sup>4</sup> for all masses
3. Central density is 5.99 3M/(4&pi;R<sup>3</sup>) for all masses

## Report (20 pts)

1. Prose is clear and at an appropriate level (**3 pts**)
2. Report acknowledges that EOS was tested (**3 pts**)
3. Discussion of convergence.  
    a. Rationale for values of &delta;, &eta;, and &xi; are provided (**5 pts**)  
    b. Scalings of central pressure, density are discussed (**3 pts**)
4. Plot shows relation with data, including error bars (**3 pts**)
5. Discussion of the realism of the white dwarf model (**3 pts**)
