using Discrete
using Base.Test

# write your own tests here
@test 1 == 1

println("Differencing tests...");
de = (n) -> 2.0^n;
for i=1:10000
  n = rand(0:100);
  
  @test de(n) == fdiff(de, n);
  @test de(n) == fdiff(de, n, order=2);
  @test de(n) == fdiff(de, n, order=3);
  @test de(n) == fdiff(de, n, order=4);
  @test de(n) == fdiff(de, n, order=rand(1:10));
  @test de(n) == bdiff(de, n);
  @test de(n) == bdiff(de, n, order=2);
  @test de(n) == bdiff(de, n, order=3);
  @test de(n) == bdiff(de, n, order=4);
  @test de(n) == bdiff(de, n, order=rand(1:10));
end

f(n)   = n^2;
df(n)  = 2*n + 1;
dbf(n) = 2*n - 1;
d2f(n) = 2;
for i=1:10000
  n = rand(0:100);
  
  @test df(n)  == fdiff(f, n);
  @test d2f(n) == fdiff(f, n, order=2);
  @test dbf(n) == bdiff(f, n);
  @test d2f(n) == bdiff(f, n, order=2);
end
println("Tests passed.");

println("Integration tests...");
for n = 4:2:8
  @test abs(simpsons(x -> 1/x, 1, e, n) - 1.0) < 10.0^(-n/2);
  @test abs(trapezoid(x -> 1/x, 1, e, n) - 1.0) < 2 * 10.0^(-n/4);

  @test abs(simpsons(x -> sin(x), -π, π, n) - 0.0) < 10.0^(-n/2);
  @test abs(trapezoid(x -> sin(x), -π, π, n) - 0.0) < 2 * 10.0^(-n/4);
  @test abs(simpsons(x -> cos(x), -π, π, n) - 0.0) < 10.0^(-n/2);
  @test abs(trapezoid(x -> cos(x), -π, π, n) - 0.0) < 2 * 10.0^(-n/4);

  @test abs(simpsons(x -> x, -1, 1, n) - 0.0) < 10.0^(-n/2);
  @test abs(trapezoid(x -> x, -1, 1, n) - 0.0) < 3 * 10.0^(-n/4);
  @test abs(simpsons(x -> x^2, -1, 1, n) - 2/3) < 10.0^(-n/2);
  @test abs(trapezoid(x -> x^2, -1, 1, n) - 2/3) < 3 * 10.0^(-n/4);
end
for n = 10:2:20
  @test abs(simpsons(x -> 1/x, 1, e, n) - 1.0) < 5e-5;
  @test abs(trapezoid(x -> 1/x, 1, e, n) - 1.0) < 5e-2;

  @test abs(simpsons(x -> sin(x), -π, π, n) - 0.0) < 5e-5;
  @test abs(trapezoid(x -> sin(x), -π, π, n) - 0.0) < 5e-2;
  @test abs(simpsons(x -> cos(x), -π, π, n) - 0.0) < 5e-5;
  @test abs(trapezoid(x -> cos(x), -π, π, n) - 0.0) < 5e-2;

  @test abs(simpsons(x -> x, -1, 1, n) - 0.0) < 5e-5;
  @test abs(trapezoid(x -> x, -1, 1, n) - 0.0) < 5e-2;
  @test abs(simpsons(x -> x^2, -1, 1, n) - 2/3) < 5e-5;
  @test abs(trapezoid(x -> x^2, -1, 1, n) - 2/3) < 5e-2;
end
println("Tests passed.");
