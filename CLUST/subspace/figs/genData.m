function GenData 
%
% GeneData.m
% Generate subspace clusters
%
% Initialize parameters
n = 10000;
k = 5;
d = 100;
Foutlier = 0;
mu = 16;
r = 2;
s = 2;

Anchor = rand(k, d)*100;
NumDim = random ('poiss', mu, 1 , k);
Dim = zeros (k , d);
index =GenRand (NumDim( 1 ) , d);
Dim(1 , find ( index ==1))=1;
for i =2: k
	dprev = min (NumDim( i -1) , floor (NumDim( i ) / 2 ) );
	down = NumDim( i ) - dprev;
	dprevindex = find ( GenRand ( dprev ,NumDim( i -1))==1);
	dprevdim = find (Dim( i-1, : ) == 1 );
	Dim( i , dprevdim ( dprevindex ) ) = 1;
	dotherdim = find (Dim( i , : ) == 0 );
	ddim = find ( GenRand ( down , d-dprev ) == 1 );
	Dim( i , dotherdim ( ddim ) ) = 1;
end

NumPoints = zeros ( 1 , k );
Nc = floor ( n*(1- Foutlier ) );
nvR = random ( 'exp' ,1 ,1 , k );
for i =1:k-1
	NumPoints ( i ) = floor(Nc*nvR ( i ) / sum( nvR ) );
end
NumPoints ( k ) = Nc-sum( NumPoints );

data = zeros ( n , d + 1 ) ;
nTemp = 0 ;
for i =1: k
	index = ( nTemp +1 ) : ( nTemp+NumPoints ( i ) ) ;
	data ( index , 1 ) = i ;
	for j =1: d
		if Dim( i , j )==0
			data ( index , j +1) = rand ( NumPoints ( i ) , 1 )* 100 ;
		else
			sij = rand *( s-1)+1;
			vRec = random ( 'norm',	Anchor ( i , j ) , sij*r , NumPoints ( i ) , 1 ) ;
			data ( index , j +1) = vRec ;
		end
	end
	nTemp = nTemp + NumPoints ( i ) ;
end
data ( ( nTemp + 1 ) : n , 1 ) = k+1;
data ( ( nTemp + 1 ) : n , 2 : ( d+1))= rand ( n-Nc , d )*100;

data = data ( randperm( n ) , : ) ;

% Pr i n t the data s e t i n f o r m a t i o n to a t e x t f i l e
fp = fopen ( '10000data100c.txt' , 'w' ) ;
fprintf( fp, 'Number of points %d, outliers %d' , n , n-Nc ) ;
for i =1: k
	fprintf ( fp , '\n Cluster %d (%d ) \n' , i , NumPoints ( i ) ) ;
	fprintf ( fp , 'Dimensions : ' ) ;
	for j =1: d
		if Dim( i , j )==1
			fprintf( fp , '%d , ' , j ) ;
		end
	end
end
fclose( fp ) ;

% Save the data s e t to a . csv f i l e
csvwrite ('10000data100c.csv' , data ) ;

%
% Generate nNum nonrepeating integers from 1 , 2 , \ l d o t s , nLen
%
function Out=GenRand (nNum, nLen )
select = zeros ( 1 , nLen ) ;
select ( floor ( rand*nLen +1)) = 1;
for i =2:nNum
	nonselect = find ( select ==0);
	nTag = nonselect ( floor ( rand*( nLen-i + 1 ) + 1 ) ) ;
	select ( nTag ) = 1 ;
end
Out = select ;

