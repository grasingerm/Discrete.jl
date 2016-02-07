module Discrete

export ConstSequence, Sequence, DFun;
export fdiff, bdiff, cdiff;
export euler, runge_kutta;

import FastAnonymous;

immutable ConstSequence
  i_start::Int;
  values::Vector{Float64};

  Sequence(i_start::Int, values::Vector{Float64}) = new(i_start, values);
end

function call(s::ConstSequence, n::Int)
  if n < s.i_start || n > s.i_start + length(s.values)
    error("Out of domain on which sequence is defined");
  end
  
  return s.values[s.i_start + n];
end 

function call(s::ConstSequence, nreal::Real)
  const n = convert(Int, nreal);
  if n < s.i_start || n > s.i_start + length(s.values)
    error("Out of domain on which sequence is defined");
  end
  
  return s.values[s.i_start + n];
end

typealias DFun     Union{Function, FastAnonymous.Fun};
typealias Sequence Union{DFun, ConstSequence};

#! Compute the forward difference of a sequence
#!
#! \param     seq     Sequence
#! \param     n       Index at which to take compute difference
#! \param     order   Order of difference
#! \return    Forward difference
function fdiff(seq::Sequence, n::Int; order::Int=1)
  @assert(order > 0, "unable to compute $order order forward difference");

  sum = 0.0;
  for i=0:order
    sum += (((order-i)%2==0) ? 
              binomial(order, i) * seq(n + i) : 
             -binomial(order, i) * seq(n + i));
  end

  return sum;
end

#! Compute the forward difference of a sequence
#!
#! \param     seq     Sequence
#! \param     n_range Range of index for which to compute differences
#! \param     order   Order of difference
#! \return    Forward differences
function fdiff(seq::Sequence, n_range::UnitRange{Int}; order::Int=1)
  return pmap(n -> fdiff(seq, n, order=order), n_range);
end

#! Compute the forward difference of a sequence
#!
#! \param     seq     Sequence
#! \param     n_range Range of index for which to compute differences
#! \param     order   Order of difference
#! \return    Forward differences
function fdiff(seq::Sequence, n_range::UnitRange{Float64}; order::Int=1)
  return pmap(n -> fdiff(seq, n, order=order), 
              convert(UnitRange{Int}, n_range));
end

#! Compute the backward difference of a sequence
#!
#! \param     seq     Sequence
#! \param     n       Index at which to take compute difference
#! \param     order   Order of difference
#! \return    Backward difference
function bdiff(seq::Sequence, n::Int; order::Int=1)
  @assert(order > 0, "unable to compute $order order forward difference");

  sum = 0.0;
  for i=0:order
    sum += (((order-i)%2==0) ? 
              binomial(order, i) * seq(n - i) : 
             -binomial(order, i) * seq(n - i));
  end

  return sum;
end

#! Compute the backward difference of a sequence
#!
#! \param     seq     Sequence
#! \param     n_range Range of index for which to compute differences
#! \param     order   Order of difference
#! \return    Backward differences
function bdiff(seq::Sequence, n_range::UnitRange{Int}; order::Int=1)
  return pmap(n -> bdiff(seq, n, order=order), n_range);
end

#! Compute the backward difference of a sequence
#!
#! \param     seq     Sequence
#! \param     n_range Range of index for which to compute differences
#! \param     order   Order of difference
#! \return    Backward differences
function bdiff(seq::Sequence, n_range::UnitRange{Float64}; order::Int=1)
  return pmap(n -> bdiff(seq, n, order=order), 
              convert(UnitRange{Int}, n_range));
end

#! Calculate forward difference on a continuous function
#!
#! \param   f     Continuous function
#! \param   x     Value to calculate difference at
#! \param   h     Step size
#! \param   order Order of difference
#! \retrun  Difference
function fdiff(f::DFun, x::Real, h::Real; order=1)
  sum = 0.0;
  for i=0:order
    sum += ((i % 2 == 0) ? binomial(order, i) * f(x + (n - i) * h) :
                          -binomial(order, i) * f(x + (n - i) * h));
  end
  return sum;
end

#! Calculate backward difference on a continuous function
#!
#! \param   b     Continuous function
#! \param   x     Value to calculate difference at
#! \param   h     Step size
#! \param   order Order of difference
#! \retrun  Difference
function bdiff(f::DFun, x::Real, h::Real; order=1)
  sum = 0.0;
  for i=0:order
    sum += ((i % 2 == 0) ? binomial(order, i) * f(x - i * h) :
                          -binomial(order, i) * f(x - i * h));
  end
  return sum;
end

#! Calculate central difference on a continuous function
#!
#! \param   b     Continuous function
#! \param   x     Value to calculate difference at
#! \param   h     Step size
#! \param   order Order of difference
#! \retrun  Difference
function cdiff(f::DFun, x::Real, h::Real; order=1)
  sum = 0.0;
  for i=0:order
    sum += ((i % 2 == 0) ? binomial(order, i) * f(x + (n/2 - i) * h) :
                          -binomial(order, i) * f(x + (n/2 - i) * h));
  end
  return sum;
end

#! Approximate the solution to an ODE using Euler's method
#!
#! \param   df    Differential equation, df = df(x,t)
#! \param   x0    Initial value of x
#! \param   ts    Range of time values     
function euler(df::DFun, x0::Real, ts::Range)
  const h = step(ts);
  x = x0;

  for t in ts
    x += h * df(x, t); 
  end

  return x;
end

#! Approximate the solution to an ODE using Runge-Kutta method
#!
#! \param   df    Differential equation, df = df(x,t)
#! \param   x0    Initial value of x
#! \param   ts    Range of time values     
function runge_kutta(df::DFun, x0::Real, ts::Range)
  const h = step(ts);
  x = x0;

  for t in ts
    k1    =   h * df(x, t);
    k2    =   h * df(x + 0.5*k1, t + 0.5*h);
    k3    =   h * df(x + 0.5*k2, t + 0.5*h);
    k4    =   h * df(x + k3, t + h);

    x    +=   1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);
  end

  return x;
end

end # end module
