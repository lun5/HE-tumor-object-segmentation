%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File: svdTutorial.m
%%%  Original Lisp Code written by D.J. Chichilnisky, summer 1992 
%%%  Converted to Matlab by D.J. Heeger, summer 1996

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This tutorial is a review of some basic concepts in linear algebra,
% concentrating on the singular value decomposition.  For a detailed
% introduction, consult a linear algebra text.  Linear Algebra and its
% Applications by Gilbert Strang is excellent.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SVD (Singular Value Decomposition):

% The SVD decomposes a matrix into the product of the three
% components:
%     A = U S V^t
% where ^t means transpose.  A is the original NxM matrix, U is
% an NxN orthonormal matrix, V is an MxM orthonormal matrix,
% and S is an NxM matrix with non-zero elements only along
% the diagonal.

% Let's try it on a randomly filled 10x5 matrix:
A = rand([10,5]);
[U,S,V] = svd(A);

% Look at S:
S

% Check that the decomposition worked:
newA = U * S * V';
mean2((newA - A).^2)

% Check that U and V are orthonormal matrices.  All of these
% should be identity matrices:
id=U'*U;
mean2((id-eye(size(id))).^2)
id=U*U';
mean2((id-eye(size(id))).^2)
id=V'*V;
mean2((id-eye(size(id))).^2)
id=V*V';
mean2((id-eye(size(id))).^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FOUR FUNDAMENTAL SUBSPACES:

% Consider a matrix A as a linear transform, y=Ax, that
% transforms N-dimension vectors, x, into M-dimensional vectors,
% y.

% If the matrix A is singular then there is a subspace of N-space
% that is mapped to zero by A, (i.e., a set of vectors x such
% that Ax=0).  This is called the "row-nullspace" of A, since
% vectors in this space are "nulled" by the row-vectors of A.

% There is a subspace of M-space that can be REACHED by A (i.e.,
% for y in this space, there exists a v in N-space such that
% Av=y).  This is called the "column-space" of A.  The
% dimensionality of the column-space is called the "rank" of A.

% The SVD explicitly constructs orthonormal bases for the
% row-nullspace and column-space of A.  The columns of U, whose
% same-numbered elements in S are non-zero, are an orthonormal
% set of basis vectors that span the column-space of A.  The
% remaining colums of U span the row-nullspace of A^t (also
% called the column-nullspace of A).

% The columns of V (rows of V^t), whose same-numbered elements in
% S are zero, are an orthonormal set of vectors that span the
% row-nullspace of A.  The remaining columns of V span the
% column-space of A^t (also called the row-space of A).

% Matlab provides functions "orth" and "null" to 

A=[1 0 0;
   0 1 0]

% Column space:
orth(A)

% Row space:
orth(A')

% Row nullspace:
null(A)

% Column nullspace (for our matrix there is no column nullspace):
null(A')


% Here's another example in which one of the columns is a linear combination
% of the others.
A =[1 -1 0; 
    0 1 -1; 
    1 0 -1]

% It is easy to see that A has rank 2 by noting that you can get
% the third column of A by summing the first 2 columns and then
% multiplying by -1.  In fact, the vector (1,1,-1) is in the
% column nullspace of A.  So is any scalar multiple of (1,1,-1).
% It is called the column nullspace because it takes the columns
% to zero: (1,1,-1) A = 0.
colnull=null(A')
colnull' * A

% The vector (1,1,1) is in the row nullspace of A.  So is any
% scalar multiple of (1,1,1).  It's called the row nullspace
% because it takes the rows to zero: A (1,1,1)^t = 0.
rownull=null(A)
A * rownull

% And the other two spaces:
orth(A)
orth(A')

% Row space and row nullspace are orthogonal:
orth(A')' * null(A)

% Col space and col nullspace are orthogonal:
orth(A)' * null(A')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SOLVING Ax=b:

A =[1 -1 0; 
    0 1 -1; 
    1 0 -1]

% The product Ax is always a combination of the columns of A:

%         (1 -1  0) (x0)      (1)      (-1)      ( 0)
%    Ax = (0  1 -1) (x1) = x0 (0) + x1 ( 1) + x2 (-1)
%         (1  0 -1) (x2)      (1)      ( 0)      (-1)

% To solve Ax=b is to find a combination of the columns that
% gives b.  We consider all possible combinations Ax, coming from
% all choices of x.  Those products form the column space of A.
% In the example, the columns lie in 3-dimensional space, but
% their combinations fill out a plane (the matrix has rank 2).
% The plane goes through the origin since one of the combinations
% has weights x0=x1=x2=0.  Some vectors b do not lie on the
% plane, so for them Ax=b can not be solved exactly.  The system
% Ax=b has an exact solution only when the right side b is in the
% column space of A.

% The vector b=(2,3,4) is not in the column space of A so Ax=b
% has no solution.  The vector b=(2,3,5) does lie in the plane
% (spanned by the columns of A) so there is a solution,
% x=(5,3,0).
A * [5 3 0]'

% However, there are other solutions as well.  x=(6,4,1) will work:
A * [6 4 1]'

% In fact, we can take x=(5,3,0) and add any scalar multiple of
% y=(1,1,1) since (1,1,1) is in the null space.  We can write it
% this way:
%    A (x + c y) = Ax + A (cy) = Ax + c Ay = Ax + c 0 = Ax

% If an NxN matrix A has linearly independent columns then
%  (1) the row nullspace contains only the point 0
%  (2) the solution to Ax=b (if there is one) is unique
%  (3) the rank of A is N
% In general any two solutions to Ax=b differ by a vector in the
% row nullspace.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% REGRESSION:

% Regression: When b is not in the column space of A, we can
% still find a vector x that comes the closest to solving
% Ax=b.

% We motivate this by thinking about fitting empirical data.  We
% change notation a bit (hopefully without too much confusion)
% and try to solve Ap=y  instead of Ax=b.

% Let's simulate collecting some data in a certain experiment.
% For each of a number of test conditions, x, we measure an
% outcome, y.  And let's assume that the relationship between y's
% and x's is given by a polynomial:
%       y = 3 x^2 + x - 4.
% It should be obvious that this was an arbitrary choice.  Let's
% use this formula to simulate some data.  First, we make a
% vector of a bunch of x  values (evenly spaced between 0 and 1).
% Then we compute y for each of those x values:
x=[0:0.1:1]';
y = 3*x.^2 + x -4;

% Look at the data:
plot(x,y,'o')

% Now let's pretend we don't know the exact relationship between
% y and x.  In particular, we have the data (xvals and yvals) and
% we know that y is a second-order polynomial of x:
%    p0 + p1 x + p2 X^2
% but we don't know the parameter values p=(p0,p1,p2)

% How might we solve for those parameter values?  We build a
% matrix A whose columns depend on the x values.  In particular,
% the first column of A has a bunch of x^2 values, the second
% column of A has a bunch of x values and the third column has a
% bunch of 1's.

N=length(x);
A=zeros(N,3);
for index=1:N
  xi=x(index);
  A(index,1)=xi^2;
  A(index,2)=xi;
  A(index,3)=1;
end

% Look at the A matrix we've just constructed:

A

% Then we solve: y = A p, where we know y and we just constructed
% A.  If A were a square (full rank) matrix, this would be easy;
% we'd just invert A to get: p = inv(A) y.  In practice, it is
% seldom necessary to form the explicit inverse of a matrix.
% Matlab provides the matrix division operator, p = A\b.  This
% produces the solution using Gaussian elimination, without
% forming the inverse.

p=A\y

% This worked perfectly because there was no noise in the data.
% Let's do it again with (simulated) noisy data.

noisyY = y + 0.1*randn(size(y));
plot(x,noisyY,'o')
p=A\noisyY

% With more data points it will work better:

x=[0:5e-3:1]';
y = 3*x.^2 + x -4;
noisyY = y + 0.1*randn(size(y));
N=length(x);
A=zeros(N,3);
for index=1:N
  xi=x(index);
  A(index,1)=xi^2;
  A(index,2)=xi;
  A(index,3)=1;
end
p=A\noisyY

% The x and y values are vectors; according to our (second-order
% polynomial) model, each element of y is a quadratic function of
% the corresponding element of xvals.  A is an Nx3 matrix, and
% the col-space of A is a 3-dimensional subspace of N space.  In
% particular, the columns of A are a basis for all possible
% second-order polynomials.  The product, Ap, is a particular
% linear combination of those basis vectors, hence, it is a
% particular second order polynomial.  If y=Ap (for some/any p),
% then we say that y is in the column space of A.  That's what we
% had in the noiseless case.  The y vector was exactly equal to
% Ap for the right choice of p.  After adding noise, however, y
% was no longer in the column space of A.  There was no choice of
% p such that y=Ap exactly.  The regression solution found a p
% such that Ap came as close as possible to y.  In other words,
% it found y-est = A p-est, where the distance between noisyY and
% y-est was minimized given that y-est was forced to be in the
% column space of A.

% All of this generalizes to many other situations, e.g., to
% models other than second-order polynials).  In particular,
% standard linear regression is a special case in which A only
% has 2 columns (x's and 1's).  We can do it for other hairier
% (non-polynomial) functions too.  Try simulating noisy data for:
% y = 10 log(x) + 20.  Plot the data.  Then estimate the
% parameters (10,20) from the noisy data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TOTAL LEAST SQUARES: 

% The regression operation we've been discussing produces a
% particular kind of "best" solution.  Let's reconsider the
% simple case of fitting a line through the origin to some (x,y)
% data points.  The regression solution finds a vector s that
% minimizes the function:
%
%      E(s) = | xs - y |^2,
%
% where s is the slope of the fitted line.  This error function
% penalizes VERTICAL displacements of the data points from the
% line.  But there are other ways to specify error functions.
% What if we believe that there is noise present in BOTH the x
% and y measurements?  We might then desire a penalty for the
% perpendicular distance of the data points from the line.  Such
% a function can easily be written as:
%
%      E'(u) = | M . u |^2 / | u |^2,
%
% where u is a 2-vector, and M is a 2-column matrix with columns
% formed from the data vectors x and y.  Since we are dividing by
% the squared norm of u, scaling u does not change the value of
% the function.  The vector that minimizes this function is
% typically scaled to be a unit vector.  It points in a direction
% perpendicular to the best-fitting line.

% This solution may be found by computing the SVD of M, and
% choosing the column of V (a 2-vector) corresponding to the
% largest singular value.  Let's try this on a data set with noise
% added to both x and y:

x=[-1:0.05:1]';
y=2*x;

noisyX = x + 0.2*randn(size(x));
noisyY = y + 0.2*randn(size(y));

plot(noisyX,noisyY,'o')

% Standard (y-distance) regression:
noisyX\noisyY

% TLS regression: (most of the time, but not always, this gives a
% better estimate)
M=[noisyX';noisyY']';
[U,S,V]=svd(M);
V(2,1)/V(1,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% COVARIANCE

% This section of the tutorial deals with multi-dimensional data
% sets.  Each data point represents one test condition, e.g., a
% data point might be the height and weight of a person
% represented as a vector (height,weight).  Or each data point
% might be the intensity values of an image represented as a very
% long vector (p1,p2,...,pN), where p_i is the intensity at the
% i_th pixel.

% The function multiRandn generates multi-dimensional random
% vector, given a particular mean and covariance matrix.  The
% following code uses multiRandn to generate a simulated data
% set, with M data points.

M=100;
ActualMean=[10 2]';
ActualCov=[1 0.8
           0.8 1];
	
data=zeros(length(ActualMean),M); 
for index=1:M
  data(:,index)=multiRandn(ActualMean,ActualCov);
end

plot(data(1,:),data(2,:),'o')

% Data is an 2xM array.  The cols of data are the data points
% that we plotted in the scatter plot.  The xvalues are in the
% first row of data, and the yvalues are in the second row.  The
% correlation coefficient in this example is 0.8, so the scatter
% plot slopes up and to the right at a 45 degree angle.

% Next let's compute the mean of the data set.  It will be close
% (but not exactly equal) to the mean that we specified at the
% outset.

EstMean=mean(data')'

% We compute the covariance of the data in two steps.  First we
% subtract EstMean from each of the columns to give us a new
% matrix D.  Then we compute D D^t, and divided by the number of
% data points, M, minus 1.

D=data-EstMean*ones(1,M);
EstCov=D*D'/(M-1)

% The built-in matlab function "cov" does the same thing (but you
% have to have the data in rows instead of columns):

cov(data')

% Obviously, the empirical covariance would be closer if we used
% a bigger data set (larger M).  Try it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PRINCIPAL COMPONENTS ANALYSIS

% Compute the svd of D.  The cols of U span the col-space of D.
% These are also the eigenvectors of the covariance, D D^t.

[U,S,V]=svd(D);
e1=U(:,1)
e2=U(:,2)

% Show that these are eigenvectors of D D^t.  For a vector v to
% be an eigenvector of a matrix A means that: Av=lv for some
% scalar l.  In other words, passing v through the matrix only
% affects its length, not its direction.

DDt = D*D';
DDt_e1 = (DDt * e1) ./ e1
DDt_e2 = (DDt * e2) ./ e2

% The eigenvector corresponding to the largest eigenvalue is
% called the first principal component.  The svd returns the
% eigenvectors in order, so e1 is the eigenvector with the
% largest eigenvalue.  The first principal component is the unit
% vector with the largest projection onto the data set.  The e1
% that you computed should be very close to (.707,.707): Note
% that you might also have gotten (-.707,-.707).  That's just as
% good.  With either sign, it points in the elongated direction
% of the scatter plot, hence it has the largest projection onto
% the data set.  The other eigenvector, e2, points in the
% perpendicular direction.  The e2 that you computed should be
% very close to (-.707,.707) or (.707,-.707):

% Again, the first principal component, e1, is the unit vector
% with the largest projection onto the data set.  The projection
% of e1 onto the data is:

norm(e1'*D)

% It's helpful to draw the shape of the matrices on a piece of
% paper: D is a (2xM) matrix, e1 is (2x1) column vector, so e1'*D
% is a (1xM) row vector .  Each element of this row vector is the
% length of the projection of a single data point onto the unit
% vector e1.  The vector-length of this row vector gives the sum
% of the projection lengths.

% The projection onto the other eigenvector is smaller:

norm(e2'*D)

% The point is that you know a lot about a data point by knowing
% only the projection of that data point onto e1.  Let's compute
% the vector-distance between a data point and its projection
% onto e1:

datapoint = D(:,1)
distance = norm(datapoint - e1*(e1'*datapoint))

% This distance is pretty small because you know a lot just by
% knowing the projection onto e1.  And now, the average of the
% vector-distances (for all the data points):

perps=D-(e1*(e1'*D));
mean(sqrt(sum(perps.^2)))

% Compare this number to what you would get by using some other
% arbitrary unit vector.  Evaluate the following a bunch of
% times.  It picks a random unit vector, prints that vector, then
% projects the data set onto that vector.  The numbers you get
% will always be bigger than or equal to what you got using e1.
% When the random unit vector is close to e1, you get a smaller
% number.  When it is far from e1, you get a larger number.

foo=rand(2,1);
e=foo/norm(foo)
perps=D-(e*(e'*D));
mean(sqrt(sum(perps.^2)))

% Try this a bunch of times to see what you get for a bunch of
% different random unit vectors.  Each data point is a vector of
% 2 numbers.  If you had to summarize each data point with one
% number, what number would you choose?

% Now do it all again for a higher dimensional data set.  Using
% this covariance matrix, what do you expect the principal
% components to be?

M=100;
ActualMean=[0 0 0]';
ActualCov=[4 0 0
           0 2 0
           0 0 1];

% Generate the data:

data=zeros(length(ActualMean),M); 
for index=1:M
  data(:,index)=multiRandn(ActualMean,ActualCov);
end

% Estimate mean and covariance:

EstMean=mean(data')'
EstCov=cov(data')

% Compute the principal components.

D=data-EstMean*ones(1,M);
[U,S,V]=svd(D);
e1=U(:,1)
e2=U(:,2)
e3=U(:,3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Local Variables:
%%% buffer-read-only: t 
%%% End:
