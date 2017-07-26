module DirectStepIntegration
using Dierckx
using DifferentialEquations
using Plots

export tkris4, tkris5, razavi

function cvs(omgd,t)
	c=cos(omgd*t);
	s=sin(omgd*t);
	return c,s	
end

function ltd_dqc(tz,p)
# Cho ra bieu dien luc tac dung dung quy cach, phu hop voi khoang thoi gian ma ket cau duoc phan tich. Vi du, luc tac dung luc dau chi duoc cho den 5 giay, nhung lai muon xet he o thoi diem 6 giay, khi do tu giay thu 5 tro di thi phai dat gia tri cua luc bang 0.
	h=tz[1];
	nstep=Int(tz[2]);
	hns=h*nstep;
	siz=size(p[2],1);
	luc=p[2];
	if size(p[1],1)==1
		ti=linspace(0,p[1]*(siz-1),siz);
	else
		ti=p[1];
	end
	tp=maximum(ti);

	if hns>tp
		noas=Int(ceil((hns-tp)/h));
		luc=[luc;zeros(noas,size(p[2],2))];
		ti=[ti;[hns-(i-1)*h for i=noas:-1:1]]
	else
		ind=searchsorted(ti,hns).stop;
		if hns>ti[ind]
			luc[ind+1,:]=luc[ind,:]+(luc[ind+1,:]-luc[ind,:])*(hns-ti[ind])/(ti[ind+1]-ti[ind]);
			ti=[ti[1:ind];hns];
			luc=luc[1:ind+1,:];
			ti=ti[1:ind+1]
		else
			luc=luc[1:ind,:];
			ti=ti[1:ind];
		end
	end
	return luc,ti
end

function tkris4(dth=([2 0;0 1],[0 0; 0 0],[96 -32;-32 32]),tz=(0.1,20),ic=([0;0],[0;0]),p=(0.1,ones(15,1)*[0 100]));
# Giai phuong trinh dao dong HE TUYEN TINH theo phuong phap Tkris
# dth: Dac trung cua he dth=(m,c,k) trong do m la ma tran khoi luong, c la ma tran he so can nhot va k la ma tran do cung
# tz: vecto co cac thong tin ve thoi gian, nhu sau tz=(h,nstep) trong do h la do lon buoc chia theo thoi gian va nstep la so buoc chia theo thoi gian.
# ic: Dieu kien ban dau ic=(cvbd,vtbd) trong do cvbd va vtbd tuong ung la vecto chuyen vi ban dau va vecto van toc ban dau
# p: Vecto tai trong p=([t_1; t_2; ...; t_n],[p(t_1);p(t_2);p(t_3);...;p(t_n)]) trong do p(t_i) la gia tri cua tai trong o thoi diem t_i. Gia tri t_1 bao gio cung bang 0. Khi cac thoi diem nay cach deu nhau, chi can vao 1 gia tri la delta_t o vecto dau tien (vi du bang 0.02).

	function xdhs_rt(ced,pse);
		hsc,hse,hsd=ced;
		ps,pe=pse;		
		rn=-((8h1m2+2(cm+3mc)*h2+2(km+6mk+4c2)*h3/5+(2ck+kc)*h4/3+h5k2/7)*hsc+(4h1mc+(3mk+c2)*h2+(4ck+kc)*h3/5+h4k2/6)*hsd+(4h1mk+h2ck+h3k2/5)*hse-(h1m1+h2c1/5+h3k1/30)*ps-(3h1m1+4h2c1/5+h3k1/6)*pe);
		sn=-((6m2+2(cm+2mc)*h+(km+3mk+3c2)*h2/2+(3ck+2kc)*h3/5+h4k2/6)*hsc+(3mc+(2mk+c2)*h+(3ck+kc)*h2/4+h3k2/5)*hsd+(3mk+h1ck+h2k2/4)*hse-(m+h1c1/4+h2k1/20)*ps-(2m+3h1c1/4+h2k1/5)*pe);	
		b=[rn;sn]
		aob=a\b;
		return aob[1:siz2], aob[siz2+1:2*siz2]
	end
	
	m,c,k=dth;
	h=tz[1];
	nstep=Int(tz[2]);
	hns=h*nstep;

	siz2=size(p[2],2);
	pt,t=ltd_dqc(tz,p);
	xxtt=[Spline1D(t,pt[:,i]) for i=1:siz2];

	u=zeros(nstep+1,siz2);
	ud=zeros(nstep+1,siz2);
	udd=zeros(nstep+1,siz2);
	
	u[1,:],ud[1,:]=ic;

	c2=c^2;
	c3=c^3;
	ck=c*k;
	cm=c*m;

	k2=k^2;	
	kc=k*c;
	km=k*m;
	
	m2=m^2;
	mc=m*c;
	mk=m*k;

	h2=h^2;
	h3=h^3;
	h4=h^4;
	h5=h^5;
	h6=h^6;
	h7=h^7;
	
	h1c1=h*c;
	h1ck=h*ck;
	h1m1=h*m;
	h1mc=h*mc;
	h1mk=h*mk;
	h1m2=h*m2;
	h2c1=h2*c;
	h2ck=h2*ck;
	h2k1=h2*k;
	h2k2=h2*k2;
	h2m2=h2*m2;
	h3k1=h3*k;
	h3k2=h3*k2;
	h3m2=h3*m2;
	h4k2=h4*k2;
	h5k2=h5*k2;
	h6k2=h6*k2;
	h7k2=h7*k2;
	
	q11=144h3m2/5+8(cm+mc)*h4+4(3km+3mk+4c2)*h5/7+(ck+kc)*h6/2+h7k2/9;
	
	q12=18h2m2+12(2cm+3mc)*h3/5+(km+2mk+2c2)*h4+(4ck+3kc)*h5/7+h6k2/8;
	q21=q12';
	
	q22=12h1m2+9(cm+mc)*h2/2+3(2km+2mk+3c2)*h3/5+(ck+kc)*h4/2+h5k2/7;
		
	a=[q11 q12; q21 q22];
	
    for i=1:nstep
    	ps=[xxtt[j](h*(i-1)) for j=1:siz2];
    	pe=[xxtt[j](h*i) for j=1:siz2];

        hse=u[i,:];
        hsd=ud[i,:];
        hsc=(m\(ps-c*hsd-k*hse))/2;
        ced=(hsc,hse,hsd);
        pse=(ps,pe);
        hsa,hsb=xdhs_rt(ced,pse);
        u[i+1,:]=hsa*h4+hsb*h3+hsc*h2+hsd*h+hse;
        ud[i+1,:]=4hsa*h3+3hsb*h2+2hsc*h+hsd;
        udd[i+1,:]=m\(pe-c*ud[i+1,:]-k*u[i+1,:]);
    end
    return u,ud,udd
end



function tkris5(dth=([2 0;0 1],[0 0; 0 0],[96 -32;-32 32]),tz=(0.1,20),ic=([0;0],[0;0]),p=(0.1,ones(15,1)*[0 100]))
# Giai phuong trinh dao dong HE TUYEN TINH theo phuong phap Tkris
# dth: Dac trung cua he dth=(m,c,k) trong do m la ma tran khoi luong, c la ma tran he so can nhot va k la ma tran do cung
# tz: vecto co cac thong tin ve thoi gian, nhu sau tz=(h,nstep) trong do h la do lon buoc chia theo thoi gian va nstep la so buoc chia theo thoi gian.
# ic: Dieu kien ban dau ic=(cvbd,vtbd) trong do cvbd va vtbd tuong ung la vecto chuyen vi ban dau va vecto van toc ban dau
# p: Vecto tai trong p=([t_1; t_2; ...; t_n],[p(t_1);p(t_2);p(t_3);...;p(t_n)]) trong do p(t_i) la gia tri cua tai trong o thoi diem t_i. Gia tri t_1 bao gio cung bang 0. Khi cac thoi diem nay cach deu nhau, chi can vao 1 gia tri la delta_t o vecto dau tien (vi du bang 0.02).

	
	function mt4(fed,pse);
		hsf,hse,hsd=fed;
		pe=pse[2];
		return mt1\(pe-k*(hsf+hse*h+hsd*h2)-c*(2hsd*h+hse)-m*(2hsd))/h;
	end;

	function mt11(fed,pse);
		hsf,hse,hsd=fed;
		return k*hsd+3c*mt4(fed,pse);
	end;

	function mt12(fed,pse);
		hsf,hse,hsd=fed;
		ps,pe=pse;
		return (k*hse+2c*hsd+6m*mt4(fed,pse)+(ps-pe)/h);
	end;

	function mt13(fed,pse);
		hsf,hse,hsd=fed;
		ps=pse[1];
		return k*hsf+c*hse+2m*hsd-ps;
	end;

	function xdhs_tk(fed,pse);
		
		sh1=(k2*mt4(fed,pse))*h6/9;
		sh2=(k'*mt11(fed,pse)+5c'*k*mt4(fed,pse))*h5/8;
		sh3=(k'*mt12(fed,pse)+5c'*mt11(fed,pse)+mt5'*k*mt4(fed,pse))*h4/7;
		sh4=(k'*mt13(fed,pse)+5c'*mt12(fed,pse)+mt5'*mt11(fed,pse)+mt6'*k*mt4(fed,pse))*h3/6;
		sh5=(5c'*mt13(fed,pse)+mt5'*mt12(fed,pse)+mt6'*mt11(fed,pse)+mt7'*k*mt4(fed,pse))*h2/5;
		sh6=(mt5'*mt13(fed,pse)+mt6'*mt12(fed,pse)+mt7'*mt11(fed,pse))*h/4;
		sh7=(mt6'*mt13(fed,pse)+mt7'*mt12(fed,pse))/3;
		sh8=(mt7'*mt13(fed,pse))/(2h);
	
		rn=-(sh1+sh2+sh3+sh4+sh5+sh6+sh7+sh8);

		sh1=(k2*mt4(fed,pse))*h5/8;
		sh2=(k'*mt11(fed,pse)+mt8'*k*mt4(fed,pse))*h4/7;
		sh3=(k'*mt12(fed,pse)+mt8'*mt11(fed,pse)+mt9'*k*mt4(fed,pse))*h3/6;
		sh4=(k'*mt13(fed,pse)+mt8'*mt12(fed,pse)+mt9'*mt11(fed,pse)+mt10'*k*mt4(fed,pse))*h2/5;
		sh5=(mt8'*mt13(fed,pse)+mt9'*mt12(fed,pse)+mt10'*mt11(fed,pse))*h/4;
		sh6=(mt9'*mt13(fed,pse)+mt10'*mt12(fed,pse))/3;
		sh7=(mt10'*mt13(fed,pse))/(2h);
	
		sn=-(sh1+sh2+sh3+sh4+sh5+sh6+sh7);
	
		b=[rn ;sn]

		aob=a\b;
		return aob[1:siz2], aob[siz2+1:2*siz2]
	
	end
	
	m,c,k=dth;
	h=tz[1];
	nstep=Int(tz[2]);
	hns=h*nstep;

	siz2=size(p[2],2);
	pt,t=ltd_dqc(tz,p);
	xxtt=[Spline1D(t,pt[:,i]) for i=1:siz2];

	u=zeros(nstep+1,siz2);
	ud=zeros(nstep+1,siz2);
	udd=zeros(nstep+1,siz2);
	
	u[1,:],ud[1,:]=ic;

	c2=c^2;
	c3=c^3;
	c4=c^4;
	c5=c^5;

	h2=h^2;
	h3=h^3;
	h4=h^4;
	h5=h^5;
	h6=h^6;
	h7=h^7;
	h8=h^8;
	h9=h^9;
	
	k2=k^2;
	k3=k^3;
	k4=k^4;
	k5=k^5;

	m2=m^2;
	m3=m^3;
	m4=m^4;

	c1k1=c*k;
	c1k2=c*k2;
	c1k3=c*k3;
	c2k1=c2*k;
	c2k2=c2*k2;
	c3k1=c3*k;
	
	kcck=k'*c+c'*k;

	mt1=6m+3h*c+h2*k;
	mt2=mt1\((-20m-5h*c-h2*k)*h2);
	mt3=mt1\((-12m-4h*c-h2*k)*h);
	mt5=(20m+k*mt2);
	mt6=(3c*mt2);
	mt7=(6m*mt2);
	mt8=(4c+k*mt3);
	mt9=(12m+3c*mt3);
	mt10=(6m*mt3);

	sh1=k2*h8/11;
	sh2=kcck*h7/2;
	sh3=(k'*mt5+mt5'*k+25c2)*h6/9;
	sh4=(k'*mt6+mt6'*k+5c'*mt5+5mt5'*c)*h5/8;
	sh5=(k'*mt7+mt7'*k+5c'*mt6+5mt6'*c+mt5'*mt5)*h4/7;
	sh6=(5c'*mt7+5mt7'*c+mt5'*mt6+mt6'*mt5)*h3/6;
	sh7=(mt6'*mt6+mt5'*mt7+mt7'*mt5)*h2/5;
	sh8=(mt6'*mt7+mt7'*mt6)*h/4;
	sh9=(mt7'*mt7)/3;

	q11=sh1+sh2+sh3+sh4+sh5+sh6+sh7+sh8+sh9;
	
	sh1=k2*h6/9;
	sh2=(k'*mt8+mt8'*k)*h5/8;
	sh3=(k'*mt9+mt9'*k+mt8'*mt8)*h4/7;
	sh4=(k'*mt10+mt10'*k+mt8'*mt9+mt9'*mt8)*h3/6;
	sh5=(mt8'*mt10+mt10'*mt8+mt9'*mt9)*h2/5;
	sh6=(mt9'*mt10+mt10'*mt9)*h/4;
	sh7=(mt10'*mt10)/3;
	
	q22=sh1+sh2+sh3+sh4+sh5+sh6+sh7;
		
	sh1=k2*h7/10;
	sh2=(k'*mt8+5c'*k)*h6/9;
	sh3=(k'*mt9+5c'*mt8+mt5'*k)*h5/8;
	sh4=(k'*mt10+5c'*mt9+mt5'*mt8+mt6'*k)*h4/7;
	sh5=(5c'*mt10+mt5'*mt9+mt6'*mt8+mt7'*k)*h3/6;
	sh6=(mt5'*mt10+mt6'*mt9+mt7'*mt8)*h2/5;
	sh7=(mt6'*mt10+mt7'*mt9)*h/4;
	sh8=(mt7'*mt10)/3;

	q12=sh1+sh2+sh3+sh4+sh5+sh6+sh7+sh8;
	q21=q12';
	
	a=[q11 q12; q21 q22];
	
    for i=1:nstep
    	ps=[xxtt[j](h*(i-1)) for j=1:siz2];
    	pe=[xxtt[j](h*i) for j=1:siz2];

        hsf=u[i,:];
        hse=ud[i,:];
        hsd=m\(ps-c*hse-k*hsf)/2;
        fed=(hsf,hse,hsd);
        pse=(ps,pe);
        hsa,hsb=xdhs_tk(fed,pse);
        hsc=mt1\((pe-m*(20hsa*h3+12hsb*h2+2hsd)-c*(5hsa*h4+4hsb*h3+2hsd*h+hse)-k*(hsa*h5+hsb*h4+hsd*h2+hse*h+hsf))/h);
        u[i+1,:]=hsa*h5+hsb*h4+hsc*h3+hsd*h2+hse*h+hsf;
        ud[i+1,:]=5hsa*h4+4hsb*h3+3hsc*h2+2hsd*h+hse;
        udd[i+1,:]=m\(pe-c*ud[i+1,:]-k*u[i+1,:]);
    end
    return u,ud,udd
end



function razavi(dth=([2 0;0 1],[0 0; 0 0],[96 -32;-32 32]),tz=(0.01,200),ic=([0;0],[0;0]),p=(0.1,ones(150,1)*[0 100]))
# Giai phuong trinh dao dong HE TUYEN TINH theo phuong phuong phap Razavi
# dth: Dac trung cua he dth=(m,c,k) trong do m la khoi luong, c la he so can nhot va k la do cung
# tz: vecto co cac thong tin ve thoi gian, nhu sau tz=[h;nstep] trong do h la do lon buoc chia theo thoi gian va nstep la so buoc chia theo thoi gian.
# ic: Dieu kien ban dau ic=(cvbd,vtbd) trong do cvbd va vtbd tuong ung la vecto chuyen vi ban dau va vecto van toc ban dau
# p: Vecto tai trong p=([t_1; t_2; ...; t_n],[p(t_1);p(t_2);p(t_3);...;p(t_n)]) trong do p(t_i) la gia tri cua tai trong o thoi diem t_i. Gia tri t_1 bao gio cung bang 0. Khi cac thoi diem nay cach deu nhau, chi can vao 1 gia tri la delta_t o vecto dau tien (vi du bang 0.02).
	function mths_ra(ced,pse);
		ps,pe=pse;
		hsc,hse,hsd=ced;
	    a=[12mh2+4ch3+kh4 6mh+3ch2+kh3; 4mh3+ch4+kh5/5  3mh2+ch3+kh4/4];
	    b=[pe-(2m+2ch+kh2)*hsc-(c+kh)*hsd-k*hse;(ps+pe)*h/2-(2mh+ch2+kh3/3)*hsc-(ch+kh2/2)*hsd-kh*hse];
	    aob=a\b;
	    return aob[1:siz2], aob[siz2+1:2*siz2]
	end
	m,c,k=dth;
	h=tz[1];
	nstep=Int(tz[2]);

	siz2=size(p[2],2);
	pt,t=ltd_dqc(tz,p);
	xxtt=[Spline1D(t,pt[:,i]) for i=1:siz2];
	
	u=zeros(nstep+1,siz2);
	ud=zeros(nstep+1,siz2);
	udd=zeros(nstep+1,siz2);
	
	u[1,:],ud[1,:]=ic;

	h2=h^2;
	h3=h^3;
	h4=h^4;
	h5=h^5;
	mh=m*h;
	mh2=m*h2;
	mh3=m*h3;
	ch=c*h;
	ch2=c*h2;
	ch3=c*h3;
	ch4=c*h4;
	kh=k*h;
	kh2=k*h2;
	kh3=k*h3;
	kh4=k*h4;
	kh5=k*h5;
	
    for i=1:nstep
       	ps=[xxtt[j](h*(i-1)) for j=1:siz2];
    	pe=[xxtt[j](h*i) for j=1:siz2];
    	
    	hsd=ud[i,:];
        hse=u[i,:];

	    hsc=(m\(ps-c*hsd-k*hse))/2;
        ced=(hsc,hse,hsd);

        pse=(ps,pe);

        hsa,hsb=mths_ra(ced,pse);
        u[i+1,:]=hsa*h4+hsb*h3+hsc*h2+hsd*h+hse;
        ud[i+1,:]=4hsa*h3+3hsb*h2+2hsc*h+hsd;
        udd[i+1,:]=m\(pe-c*ud[i+1,:]-k*u[i+1,:]);
    end
    return u,ud,udd
end

function newmark(dth=([2 0;0 1],[0 0; 0 0],[96 -32;-32 32]),tz=(0.01,200),ic=([0;0],[0;0]),p=(0.1,ones(150,1)*[0 100]),param=(1/6,1/2))

# Giai phuong trinh dao dong HE TUYEN TINH theo phuong phuong phap Newmark
# dth: Dac trung cua he dth=(m,c,k) trong do m la khoi luong, c la he so can nhot va k la do cung
# tz: vecto co cac thong tin ve thoi gian, nhu sau tz=[h;nstep] trong do h la do lon buoc chia theo thoi gian va nstep la so buoc chia theo thoi gian.
# ic: Dieu kien ban dau ic=(cvbd,vtbd) trong do cvbd va vtbd tuong ung la vecto chuyen vi ban dau va vecto van toc ban dau
# p: Vecto tai trong p=([t_1; t_2; ...; t_n],[p(t_1);p(t_2);p(t_3);...;p(t_n)]) trong do p(t_i) la gia tri cua tai trong o thoi diem t_i. Gia tri t_1 bao gio cung bang 0. Khi cac thoi diem nay cach deu nhau, chi can vao 1 gia tri la delta_t o vecto dau tien (vi du bang 0.02).

	m,c,k=dth;
	ps=p[1];
	h=tz[1];
	nstep=Int(tz[2]);
	h2=h^2;
	beta,gamma=param;
	bh=beta*h;
	bh2=beta*h2;
	hg2b=h*(gamma/(2beta)-1);
	gobh=gamma/bh;
	gob=gamma/beta;

	siz2=size(p[2],2);
	pt,t=ltd_dqc(tz,p);
	xxtt=[Spline1D(t,pt[:,i]) for i=1:siz2];

	u=zeros(nstep+1,siz2);
	ud=zeros(nstep+1,siz2);
	udd=zeros(nstep+1,siz2);
	
	u[1,:],ud[1,:]=ic;

    udd[1,:]=m\(p[2][1,:]-c*ud[1,:]-k*u[1,:]);
    k_hat=k+gobh*c+m/bh2;
    a=m/bh+gob*c;
    b=m/(2beta)+hg2b*c;
    for i=1:nstep
       	ps=[xxtt[j](h*(i-1)) for j=1:siz2];
    	pe=[xxtt[j](h*i) for j=1:siz2];

        delta_p=pe-ps;
        delta_p_hat=delta_p+a*ud[i,:]+b*udd[i,:];
        delta_u=k_hat\delta_p_hat;
        delta_udot=gobh*delta_u-gob*ud[i,:]-hg2b*udd[i,:];
        delta_uddot=delta_u/bh2-ud[i,:]/bh-udd[i,:]/(2beta);
        u[i+1,:]=u[i,:]+delta_u;
        ud[i+1,:]=ud[i,:]+delta_udot;
        udd[i+1,:]=udd[i,:]+delta_uddot;
    end
    return u,ud,udd
end

function wilson(dth=([2 0;0 1],[0 0; 0 0],[96 -32;-32 32]),tz=(0.01,200),ic=([0;0],[0;0]),p=(0.1,ones(150,1)*[0 100]),param=1.4)
# Giai phuong trinh dao dong HE TUYEN TINH theo phuong phuong phap Newmark
# dth: Dac trung cua he dth=[m;c;k] trong do m la khoi luong, c la he so can nhot va k la do cung
# tz: vecto co cac thong tin ve thoi gian, nhu sau tz=[h;nstep] trong do h la do lon buoc chia theo thoi gian va nstep la so buoc chia theo thoi gian.
# ic: Dieu kien ban dau ic=[cvbd;vtbd] trong do cvbd va vtbd tuong ung la chuyen vi ban dau va van toc ban dau
# p: Vecto tai trong p=[p(0);p(t1);p(t2);...;p(t_n)] trong do p(t_i) la gia tri cua tai trong o thoi diem t_i
# param: khai bao he so theta cua phuong phap
	m,c,k=dth;
	ps=p[1];
	h=tz[1];
	nstep=Int(tz[2]);
	h2=h^2;
	theta=param;
	thh=theta*h;
	thh2=thh^2;
    hs1=6m/thh2+3c/thh+k;

	siz2=size(p[2],2);
	pt,t=ltd_dqc(tz,p);
	xxtt=[Spline1D(t,pt[:,i]) for i=1:siz2];

	u=zeros(nstep+1,siz2);
	ud=zeros(nstep+1,siz2);
	udd=zeros(nstep+1,siz2);
	
	u[1,:],ud[1,:]=ic;

    udd[1,:]=m\(p[2][1,:]-c*ud[1,:]-k*u[1,:]);

    for i=1:nstep
       	ps=[xxtt[j](h*(i-1)) for j=1:siz2];
    	pe=[xxtt[j](h*i) for j=1:siz2];
    	pith=ps*(1-theta)+pe*theta;
        shtd1=m*(6u[i,:]/thh2+6ud[i,:]/thh+2udd[i,:])+c*(3u[i,:]/thh+2ud[i,:]+thh*udd[i,:]/2)+pith;
        uith=hs1\shtd1;
        uddith=6*(uith-u[i,:]-thh*ud[i,:]-thh2*udd[i,:]/3)/thh2;
        udd[i+1,:]=udd[i,:]+(uddith-udd[i,:])/theta;
        ud[i+1,:]=ud[i,:]+h*(udd[i,:]+udd[i+1,:])/2;
        u[i+1,:]=u[i,:]+h*ud[i,:]+h2*udd[i,:]/3+h2*udd[i+1,:]/6;
    end
    return u,ud,udd
end
end
