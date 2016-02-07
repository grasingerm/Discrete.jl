using Discrete
using Base.Test

# write your own tests here
@test 1 == 1

e = (n) -> 2.0^n;
for i=1:10000
  n = rand(0:1000);
  
  e(n) == fdiff(e, n);
  e(n) == fdiff(e, n, order=2);
  e(n) == fdiff(e, n, order=3);
  e(n) == fdiff(e, n, order=4);
  e(n) == fdiff(e, n, order=rand(1:25));
  e(n) == bdiff(e, n);
  e(n) == bdiff(e, n, order=2);
  e(n) == bdiff(e, n, order=3);
  e(n) == bdiff(e, n, order=4);
  e(n) == bdiff(e, n, order=rand(1:25));
end

f(n)   = n^2;
df(n)  = 2*n + 1;
dbf(n) = 2*n - 1;
d2f(n) = 2;
for i=1:10000
  n = rand(0:1000);
  
  df(n)  == fdiff(e, n);
  d2f(n) == fdiff(e, n, order=2);
  dbf(n) == bdiff(e, n);
  d2f(n) == bdiff(e, n, order=2);
end
