%Function for bluring the image
%This will output the image with a pixel value equal to the mean of the
%pixels in a submatrix 2w+1 where the given pixel sits in the middle. e.g.
%if w = 1 then the submatrix will be 3x3, if w=2 then it will be 5x5 e.c.t.

function output = blur(img,w)
B=double(img);
[m,n] = size(B);
k=2*w+1;
for i = 1:m %i is equal to number of horizontal pixels
    for j = 1:n %j is equal to number of veritcal pixels
        p=i-fix(k/2);
        q=i+fix(k/2);
        r=j-fix(k/2);
        s=j+fix(k/2);
        if p<1
            p=1;
        end
        if q>m
            q=m;
        end
        if r<1
            r=1;
        end
        if s>n
            s=n;
        end
        A=B([p:q],[r:s]);
        C(i,j)=mean(A(:));
    end
end
output=uint8(C);
end