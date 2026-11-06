function [bestValue,bestX,curve]=soa_paper_engine(N,M,lb,ub,D,fobj,initialKeys)
% Dhiman/Kumar 2019, equations (5)-(14), fc=2,u=v=1.
% rd/k: independent per agent-coordinate; global incumbent, box clipping.
X=lb+(ub-lb)*rand(N,D); if ~isempty(initialKeys), X(1,:)=initialKeys; end
bestValue=Inf; bestX=zeros(1,D); curve=zeros(M,1);
for generation=0:M-1
    X=max(lb,min(ub,X));
    for i=1:N
        value=fobj(X(i,:)); if value<bestValue, bestValue=value; bestX=X(i,:); end
    end
    A=2-generation*2/M;
    for i=1:N
        rd=rand(1,D); k=2*pi*rand(1,D); B=2*A^2*rd;
        Cs=A*X(i,:); Ms=B.*(bestX-X(i,:)); Ds=abs(Cs+Ms);
        radius=exp(k); xx=radius.*cos(k); yy=radius.*sin(k); zz=radius.*k;
        X(i,:)=Ds.*xx.*yy.*zz+bestX;
    end
    curve(generation+1)=bestValue;
end
end
