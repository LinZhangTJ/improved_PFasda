%transform the SHC vector of 3717*1 to CS 61*61
function cs=vector2cs(vector,N)
sz=0;
for iin=3:N+1
    num=2*(iin-1)+1;
    sc(iin,N+1-(num-1)/2:N+1+(num-1)/2)=vector(1+sz:num+sz);
    sz=sz+num;
end
cs=gmt_sc2cs(sc);