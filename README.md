# Brunner_Munzel_CPP
 C++ Statistical Tests: Brunner Munzel & Welch t


## Usage

```c++:main
/* 2 datasets */  
std::vector<double> x = { 1,2,1,1,1,1,1,1,1,1,2,4,1,1 };  
std::vector<double> y = { 3,3,4,3,1,2,3,1,1,5,4 };

/* brunner munzel test */  
test::test_base<double>& res1 = test::brunner_munzel<double>(x, y);  
printf("g: %lf, l:%lf, b:%lf\n", res1.pvalue_greater(), res1.pvalue_less(), res1.pvalue_both());

/* welch t test */  
test::test_base<double>& res2 = test::welch_t<double>(x, y);  
printf("g: %lf, l:%lf, b:%lf\n", res2.pvalue_greater(), res2.pvalue_less(), res2.pvalue_both());
```

