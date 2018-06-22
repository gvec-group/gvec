#!/bin/bash

../build/bin/gvec parameter_test.ini |tee log.run_gvec

../build/bin/gvec parameter_test.ini GVEC_TEST_State_0000_99999999.dat |tee log.restart_gvec

../build/bin/test_gvec_to_hopr GVEC_TEST_State_0000_99999999.dat |tee log.tohopr

../build/bin/test_gvec_to_gene GVEC_TEST_State_0000_99999999.dat |tee log.togene
