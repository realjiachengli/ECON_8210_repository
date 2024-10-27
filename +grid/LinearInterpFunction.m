classdef LinearInterpFunction < grid.ApproxFunction
    
    properties (SetAccess=protected)
        % inherited properties (abstract in superclass)
        SSGrid
        Nof
        Vals      
        % linear interp specific properties
        Type % scatter or equidistant
        InterpStruct
        SingletonDim
    end
       
    methods
        % constructor
        function sf=LinearInterpFunction(ssgrid,vals)
            sf.SSGrid=ssgrid;
            sf=fitTo(sf,vals);
        end
        
        function sf=set.SSGrid(sf,ssg)
            % check that grid object at initialization is TensorGrid
            if ~isa(ssg,'grid.TensorGrid')
                error('StateSpaceGrid must be a TensorGrid or a ScatterGrid');
            end
            sf.SSGrid=ssg;
        end
        

        
        % fit to new values
        function sf=fitTo(sf,vals)
            [npt,nof]=size(vals);
            if npt~=sf.SSGrid.Npt
                error('Value matrix must have dimensions (Npt x Nof)');
            end
            sf.Nof=nof;
            sf.Vals=vals;
            % reshape values for call to spapi
            dimvec=sf.SSGrid.Dimvec';
            ndim=sf.SSGrid.Ndim;
            singleton=(dimvec==1);
            if  sum(singleton)>0
                dimvec=dimvec(~singleton);
                ndim=ndim-sum(singleton);
                sf.SingletonDim=singleton;
            else
                sf.SingletonDim=[];
            end
            unigrids=sf.SSGrid.Unigrids(~singleton);
            vals=reshape(vals,[dimvec,nof]);
            vals=permute(vals,[ndim+1,1:ndim]);
            baseindex=repmat({':'},1,ndim+1);
            sf.InterpStruct=cell(nof,1);
            for f=1:nof
                index=baseindex;
                index{1}=f;
                sf.InterpStruct{f}=griddedInterpolant(unigrids,squeeze(vals(index{:})),'linear','none');
            end
            
        end
        
        % evaluation
        function vals=evaluateAt(sf,points)
            [np,ndim]=size(points);
            if ndim~=sf.SSGrid.Ndim
                error('Point matrix must have dimensions (#points x Ndim)');
            end
            % extrapolation doesn't work well, so force back into state
            % bounds
            points=real(points);
            points_corr=points;
            SBlow=ones(np,1)*sf.SSGrid.StateBounds(1,:);
            SBhi=ones(np,1)*sf.SSGrid.StateBounds(2,:);
            upvio=(points>SBhi);
            points_corr(upvio)=SBhi(upvio);
            downvio=(points<SBlow);
            points_corr(downvio)=SBlow(downvio);
            if ~isempty(sf.SingletonDim)
                ndim=size(points_corr,2);
                colind=(1:ndim);
                points_corr=points_corr(:,colind(~sf.SingletonDim));
            end           
            vals=zeros(sf.Nof,np);
            for f=1:sf.Nof
                vals(f,:)=sf.InterpStruct{f}(points_corr);
            end
        end
        
        function plot2D(cf,dispfct,dispdim,val_other,ln,fig)         
			if nargin<6
			   fig=true;
			end
		   
           totdim=cf.SSGrid.Ndim;
           dim_other=setdiff(1:totdim,dispdim);
           
           gridx=cf.SSGrid.Unigrids{dispdim}';
           np=length(gridx);
           
           plist=zeros(np,totdim);
           plist(:,dispdim)=gridx;
           for i=1:length(val_other)
               plist(:,dim_other(i))=ones(np,1)*val_other(i);
           end
           
           vlist=evaluateAt(cf,plist);
           vlist=vlist(dispfct,:);
           
           if ln==1
               vlist=exp(vlist);
           elseif ln==2
               vlist=1./(1+exp(vlist));
		   end
		   if fig
               figure;
		   end
           plot(gridx,vlist,'.-');
       end        
    
        
    end
    
    
end
