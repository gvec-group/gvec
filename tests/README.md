## Tests 

### Running these tests

After having compiled gvec inot the build folder, run here:

``` 
   ../build/bin/gvec parameter_test.ini
```
If successfull,  test the restart
``` 
   ../build/bin/gvec parameter_test_restart.ini GVEC_TEST_State_0000_99999999.dat
``` 
and also run the tests for the interfaces
``` 
   ../build/bin/test_gvec_to_hopr GVEC_TEST_State_0000_99999999.dat
   ../build/bin/test_gvec_to_gene GVEC_TEST_State_0000_99999999.dat
``` 

