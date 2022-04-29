# Contributing
The development of the analysis package is more than welcome. So feel free to contribute with your ideas

## Contributing procedure
I inspire you one to use separate branch for the development of your analysis. Thus, you can incorporate changes made to the develop branch to your analysis at any time. If you want to contribute to dev branch, create a marge request and ask for the review, so the people can see what exactly you changed and how it may affect other analysis.

## Styling convention
Please follow the simple code convention we have in the package so far:
1. leave an underscore `_` before the class variable
2. Comment your methods. The best way is to follow `///` Doxygen comment convention for the short comment and `/** */` for the long one
3. Line length < 80 is preferable

## Best practice:
1. Use numeric cast from Geom.hxx `num::cast<float>(a)`
to control the narrowing and signess. 

## Others
The separate [testbeam-plotter](https://gitlab.com/t2k-beamtest/testbeam_plotters) package is prefered for macros upload