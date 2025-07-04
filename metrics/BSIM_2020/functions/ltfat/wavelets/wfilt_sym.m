function [h,g,a,info]=wfilt_sym(N)
%WFILT_SYM Symlet filters 
%   Usage: [h,g,a]=wfilt_sym(N);    
%
%   [h,g,a]=WFILT_SYM(N) generates the "least asymmetric" Daubechies'
%   orthogonal wavelets or "symlets" with N vanishing moments and 
%   length 2N.  
%   Zeros of the trigonometrical polynomial the filters consist of in the 
%   Z-plane are selected alternatingly inside and outside the unit circle.
%
%   Remark: Filters generated by this routine differ slightly from the
%   ones in the reference (table 6.3, figure. 6.4) because of the ambiguity
%   in the algorithm.
%
%   Examples:
%   ---------
%   :
%     wfiltinfo('sym8');
%   
%   References:
%     I. Daubechies. Ten Lectures on Wavelets. Society for Industrial and
%     Applied Mathematics, Philadelphia, PA, USA, 1992.
%     
%     
%
%   Url: http://ltfat.github.io/doc/wavelets/wfilt_sym.html

% Copyright (C) 2005-2015 Peter L. Soendergaard <peter@sonderport.dk>.
% This file is part of LTFAT version 2.1.1
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Original copyright goes to:
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es

num_coefs = 2*N;
a = [2;2];
info.istight = 1;

if num_coefs==2	    % Haar filters
	[h,g,a,info]=wfilt_db(1);
	return
end

N=num_coefs/2;

poly=trigpol(N);    %Calculate trigonometric polynomial 

ceros=roots(poly);  %Calculate roots

realzeros=[];
imagzeros=[];
numrealzeros=0;
numimagzeros=0;

for i=1:2*(N-1)
	if (imag(ceros(i))==0)
		numrealzeros=numrealzeros+1;
		realzeros(numrealzeros)=ceros(i);
	else
		numimagzeros=numimagzeros+1;
		imagzeros(numimagzeros)=ceros(i);
	end
end


%% complex zeros are grouped together
i=0;
for cont=1:numimagzeros/4
	modulos(cont)=abs(imagzeros(cont+i));
	alfa(cont)=angle(imagzeros(cont+i));
	i=i+1;
end

%% Calculate phase contribution of complex and real zeros for all the
%% combination of these zeros. Each group of zeros is identified with a binary
%% number.

indice=2^(numimagzeros/4+numrealzeros/2);
fase=zeros(indice,1001);
for cont=0:indice-1,
	bin=dec2bina(cont,log2(indice));
   	for i=1:length(bin)-numrealzeros/2
		if bin(i)
			R=1/modulos(i);
		else
			R=modulos(i);
		end
		alf=alfa(i);
		fase(cont+1,:)=fase(cont+1,:)+atang1(R,alf);
	end
	ind=1;
	for i=length(bin)-numrealzeros/2+1:length(bin)
		if bin(i)
			R=realzeros(ind+1);		
			R=realzeros(ind+1);
		else
			R=realzeros(ind);
		end
		ind=ind+2;
	 	fase(cont+1,:)=fase(cont+1,:)+atang2(R);

	end	
end

%% To retain only the non linear part of the phase.

fas=linealiz(fase);

imagzeros=[];
zerosreales=[];


%% To see which phase is closer to zero we select the one with minimun variance

[maximo,pos]=min(sum(fas'.^2));  

bin=dec2bina(pos-1,log2(indice));

for i=1:length(bin)-numrealzeros/2
	if bin(i)
		z1=1/modulos(i)*exp(j*alfa(i));
	else
		z1=modulos(i)*exp(j*alfa(i));	
	end
	imagzeros=[imagzeros z1 conj(z1)];
end

ind=1;
for i=length(bin)-numrealzeros/2+1:length(bin)
	if bin(i)
		zerosreales=[zerosreales realzeros(ind+1)];
	else
		zerosreales=[zerosreales realzeros(ind)];
	end
	ind=ind+2;
end

% Construction of rh from its zeros

numrealzeros=numrealzeros/2;
numimagzeros=numimagzeros/2;

rh=[1 1];

for i=2:N
	rh=conv(rh,[1 1]);
end

for i=1:numrealzeros
	rh=conv(rh,[1 -zerosreales(i)]);
end

for i=1:2:numimagzeros
	rh=conv(rh,[1 -2*real(imagzeros(i)) abs(imagzeros(i))^2]);
end

% Once ho is factorized in its zeros, it must be normalized multiplying by "cte".

cte=sqrt(2)/sum(rh);
rh=cte*rh;
fLen = length(rh);

% Some odd values of N produce flipped filters
% Bigger N jut take forever to calculate.
if any(N==[7,9]) || ( N>=13 && rem(N,2) == 1)
    rh = rh(end:-1:1);
end

g{1} = rh;
g{2} = -(-1).^(0:fLen-1).*g{1}(end:-1:1);
Lh = numel(rh);
d = cellfun(@(gEl) -length(gEl)/2,g);
if N>2
  % Do a filter alignment according to "center of mass"
  d(1) = -find(abs(g{1})==max(abs(g{1})),1,'first')+1;
  d(2) = -find(abs(g{2})==max(abs(g{2})),1,'first')+1;
  if abs(rem(d(1)-d(2),2))==1
      % Shift d(2) just a bit
      d(2) = d(2) - 1;
  end
end

g = cellfun(@(gEl,dEl) struct('h',gEl(:),'offset',dEl),g,num2cell(d),'UniformOutput',0);
 
h = g;



function  bin=dec2bina(num,bits)

%DEC2BINA    BIN = DEC2BINA(NUM,BITS) returns a vector which contains 
%	     the decimal number NUM in binary format, with a number of 
%	     digits equal to BITS. It is an auxiliary function used by
%	     SYMLETS.

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
%                                                      
%                                                      
% Uvi_Wave is free software; you can redistribute it and/or modify it      
% under the terms of the GNU General Public License as published by the    
% Free Software Foundation; either version 2, or (at your option) any      
% later version.                                                           
%                                                                          
% Uvi_Wave is distributed in the hope that it will be useful, but WITHOUT  
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    
% for more details.                                                        
%                                                                          
% You should have received a copy of the GNU General Public License        
% along with Uvi_Wave; see the file COPYING.  If not, write to the Free    
% Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.             
%                                                                          
%       Author: Jose Martin Garcia
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------


if nargin<2
	flag=0;
else
	flag=1;
end

bin=[];
coc=num;
while coc>1
	bin=[rem(coc,2) bin];
	coc=fix(coc/2);
end
bin=[coc bin];
if flag 
 	if length(bin)<bits
		bin=[zeros(1,bits-length(bin)) bin];
	end
end

function fase=atang1(R,alfa)

%ATANG1    PHASE=ATANG1(R,ALFA) returns the phase contribution
%	   of a complex pair of zeros. Linear terms have been
%	   removed. It is an auxiliary function used by SYMLETS.

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
%                                                      
%                                                      
% Uvi_Wave is free software; you can redistribute it and/or modify it      
% under the terms of the GNU General Public License as published by the    
% Free Software Foundation; either version 2, or (at your option) any      
% later version.                                                           
%                                                                          
% Uvi_Wave is distributed in the hope that it will be useful, but WITHOUT  
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    
% for more details.                                                        
%                                                                          
% You should have received a copy of the GNU General Public License        
% along with Uvi_Wave; see the file COPYING.  If not, write to the Free    
% Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.             
%                                                                          
%       Author: Jose Martin Garcia
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------


w=[0:2*pi/1e3:2*pi];		%frequency axis

fas=atan( (1-R^2)*sin(w)./((1+R^2)*cos(w)-2*R*cos(alfa)) );

zero=acos(2*R*cos(alfa)/(1+R^2));
u1=ceil(zero*1000/(2*pi))+1;
u2=ceil((2*pi-zero)*1000/(2*pi))+1;
if (1-R^2)*sin(zero)<0
	cte=pi;
	fase=fas+w;
else
	fase=fas-w;
	cte=-pi;
end
fase(u1:1001)=fase(u1:1001)-cte;
fase(u2:1001)=fase(u2:1001)-cte;

function fase=atang2(R)

%ATANG2    PHASE=ATANG2(R) returns the phase contribution of
%	   a real zero. Linear terms have been removed. It is
%	   an auxiliary function used by SYMLETS.

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
%                                                      
%                                                      
% Uvi_Wave is free software; you can redistribute it and/or modify it      
% under the terms of the GNU General Public License as published by the    
% Free Software Foundation; either version 2, or (at your option) any      
% later version.                                                           
%                                                                          
% Uvi_Wave is distributed in the hope that it will be useful, but WITHOUT  
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    
% for more details.                                                        
%                                                                          
% You should have received a copy of the GNU General Public License        
% along with Uvi_Wave; see the file COPYING.  If not, write to the Free    
% Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.             
%                                                                          
%       Author: Jose Martin Garcia
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------


w=[0:2*pi/1e3:2*pi];	%frequency axis

fas=atan( (1+R)/(1-R)*tan(w/2) );

if R<1
	fase=fas-w;
	cte=-pi;
else
	fase=fas+w;
	cte=pi;
end;
u=ceil(pi*1000/(2*pi))+2;
fase(u:1001)=fase(u:1001)-cte;

function fase=linealiz(f)

%LINEALIZ     PHASE = LINEALIZ(F) is an auxiliary function used
%	      by SYMLETS. It eliminates the linearity of the
%	      phase vector F.


%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
%                                                      
%                                                      
% Uvi_Wave is free software; you can redistribute it and/or modify it      
% under the terms of the GNU General Public License as published by the    
% Free Software Foundation; either version 2, or (at your option) any      
% later version.                                                           
%                                                                          
% Uvi_Wave is distributed in the hope that it will be useful, but WITHOUT  
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    
% for more details.                                                        
%                                                                          
% You should have received a copy of the GNU General Public License        
% along with Uvi_Wave; see the file COPYING.  If not, write to the Free    
% Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.             
%                                                                          
%       Author: Jose Martin Garcia
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------


if abs(sum(f(1,:))) >0.0001
	w=[0:2*pi/1e3:2*pi];
	[m,n]=size(f);
	fase=zeros(m,n);
	for cont=1 : m
		if sum(f(cont,:)) >0
			fase(cont,:)=f(cont,:)-w/2;
		else
			fase(cont,:)=f(cont,:)+w/2;
		end
	end
else
	fase=f;
end 

function polinomio=trigpol(N)

coefs=zeros(N,2*N-1);
coefs(1,N)=1;

 
for i=1:N-1
	fila=[1 -2 1];
	for j=2:i
		fila=conv(fila,[1 -2 1]);
	end;
	fila=numcomb(N-1+i,i)*(-0.25)^i*fila;
	fila=[ zeros(1,(N-i-1))  fila zeros(1,(N-i-1))];
	coefs(i+1,:)=fila;
end

for i=0:(2*(N-1))
	polinomio(i+1)=0;
	for j=1:N
		polinomio(i+1)=polinomio(i+1)+coefs(j,i+1);
	end
end; 

function y=numcomb(n,k)

if n==k,
   y=1;
elseif k==0,
   y=1;
elseif k==1,
   y=n;
else 
   y=fact(n)/(fact(k)*fact(n-k));
end

function y=fact(x)

% FACT   Factorial.
%        FACT(X) is the factorial of the elements in X vector.

for j=1:length(x)
    if x(j)==0,
       y(j)=1;
    else
       y(j)=x(j)*fact(x(j)-1);
    end
end




