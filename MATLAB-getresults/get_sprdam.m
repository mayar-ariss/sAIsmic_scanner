function [ DAMS,DAMO ] = get_sprdam(SPRP,EN,DEN,file,fname,DIM_EL,nst,plot_dams,plot_damo)
%Obtain spring behaviour from num file and derive damage data for shear
%springs (DAMS) and for out-of-plane springs (DAMO)

DAMS=cell(length(file),1);
DAMO=cell(length(file),1);

for pp=1:length(file)
    % Get spring properties

	%In-plane spring
	kpl1=SPRP{pp}(11);
	fy=SPRP{pp}(12);
    fult=SPRP{pp}(18); %used for both IP and OOP spring
	Ed=SPRP{pp}(1);
	Uult=SPRP{pp}(4);
	
	%Out-of-plane spring
	py=SPRP{pp}(6);
	pr=SPRP{pp}(7);
	muy=SPRP{pp}(8);
	mur=SPRP{pp}(9);
	Kp=4*SPRP{pp}(10);
	%U=SPRP{pp}(20);
	U=4;
    
    % Get all necessary data for all the springs from each partition num
    % file
    if strcmp(DEN{pp},'all')
        DEN1=EN{pp};
    else
        DEN1=DEN{pp};
    end
    FUD = get_FUdiag25(DEN1,file(pp),fname,nst);

	M=zeros(nst,((length(DEN1)*2)+1));
	M(:,1)=[1:1:nst];
		
    for ee=1:size(DEN1,1);
		
		B=DIM_EL{pp}(ee,1);
		H=DIM_EL{pp}(ee,2);
		Uult_spring=Uult*cos(atan(H/B));
		x=3;
		y=4;
        % calculate damage of each element based on f-u hisotry
        Fd=FUD{ee}(:,1); Ud=FUD{ee}(:,2);
        Fo=FUD{ee}(:,3); Uo=FUD{ee}(:,4);
		
	if plot_dams==0
		M(:,ee*2)=Ud;
		M(:,ee*2+1)=Fd;

		m=size(DEN1,1);
		xb=[0;0;0;0;0;0;0];
		yb=[0;0;0;0;0;0;0];
        %  shear diagonal
		for i=1:(nst-4)
			if (((Fd(y)/Ud(y))/(Fd(x)/Ud(x)))>0.999)&&(((Fd(y)/Ud(y))/(Fd(x)/Ud(x)))<1.001)
				x=x+1;
				y=y+1;
			else
				break
			end
		end
		if y==4+(nst-4)
			DAMS{pp}(ee)=0;
			%hold off
		else
		Fy=abs(Fd(x)+Fd(y))/2;
		uy=abs((Ud(x)+Ud(y))/2);
		Fy_neg=-Fy;
		uy_neg=-uy;

		Kel=Fy/uy;
		Kpl=Kel*kpl1;

		Fmax=Fy/fy;
		umax=uy+(Fmax-Fy)/Kpl;
		Fmax_neg=-Fmax;
		umax_neg=-umax;
	
		Fult=0.5*Fmax;
		Fult_neg=-Fult;
	
		Kpost=-Fmax/(Uult_spring-umax);
		
		uult=abs((Fult-Fmax)/Kpost+umax);
		uult_neg=-uult;
		
		xb=[uult_neg;umax_neg;uy_neg;0;uy;umax;uult];
		yb=[Fult_neg;Fmax_neg;Fy_neg;0;Fy;Fmax;Fult];
				
		upos=max(Ud);
		uneg=min(Ud);
		u=max(abs(upos),abs(uneg));
		uy=abs(uy);

        DAMS{pp}(ee)=(u-uy)/(Uult_spring-uy);
			
		end
end 
		

if plot_damo==0
		M(:,ee*2)=Uo;
		M(:,ee*2+1)=Fo;
		x=3;
		y=4;
        % calculate damage of each element based on f-u hisotry
        Fo=FUD{ee}(:,3); Uo=FUD{ee}(:,4);
		m=size(DEN1,1);
		
		xb=[0;0;0;0;0;0;0];
		yb=[0;0;0;0;0;0;0];
        %  out-of-plane diagonal shear spring 
		for i=1:(nst-4)
			if (((Fo(y)/Uo(y))/(Fo(x)/Uo(x)))>0.999)&&(((Fo(y)/Uo(y))/(Fo(x)/Uo(x)))<1.001)
				x=x+1;
				y=y+1;
			else
				break
			end
		end
		if y==4+(nst-4)
			DAMO{pp}(ee)=0;
		else
		Fy_o=abs(Fo(x)+Fo(y))/2;
		uy_o=abs((Uo(x)+Uo(y))/2);
		Fy_neg_o=-Fy_o;
		uy_neg_o=-uy_o;

		Kel_o=Fy_o/uy_o;
			
		sdt=abs(((Fy_o/(2*B*H))-py)/muy);
		Fr_o=(pr+mur*sdt)*2*H*B;
		Fr_neg_o=-Fr_o;
		ur_o=uy_o+(Fr_o-Fy_o)/Kp;
		ur_neg_o=-ur_o;
		uult_o=U*ur_o;
		uult_neg_o=-ur_o;

		xb=[uult_neg_o;ur_neg_o;uy_neg_o;0;uy_o;ur_o;uult_o];
		yb=[Fr_neg_o;Fr_neg_o;Fy_neg_o;0;Fy_o;Fr_o;Fr_o];
						
		upos_o=max(Uo);
		uneg_o=min(Uo);
		u_o=max(abs(upos_o),abs(uneg_o));
		uy_o=abs(uy_o);
        if u_o<uy_o
            DAMO{pp}(ee)=0;
		elseif u_o>ur_o
			DAMO{pp}(ee)=1;
        else
			DAMO{pp}(ee)=(u_o-uy_o)/(ur_o-uy_o);
		end	   

    end
end
end
end

