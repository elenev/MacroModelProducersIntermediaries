classdef LinearInterpFunction < grid.ApproxFunction
    
    properties (SetAccess=protected)
        % inherited properties (abstract in superclass)
        SSGrid
        Nof
        Vals      
        % linear interp specific properties
        Type % scatter or equidistant
        InterpStruct
    end
    
    properties (Constant)
        SCATTER=10;
        EQUI=20;
    end
    
    methods
        % constructor
        function sf=LinearInterpFunction(ssgrid,vals)
            sf.SSGrid=ssgrid;
            sf=fitTo(sf,vals);
        end
        
        function sf=set.SSGrid(sf,ssg)
            % check that grid object at initialization is TensorGrid
            if isa(ssg,'grid.TensorGrid')
                sf.Type=grid.LinearInterpFunction.EQUI;
            elseif isa(ssg,'grid.ScatterGrid')
                sf.Type=grid.LinearInterpFunction.SCATTER;
            else
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
            % distinguish cases
            if sf.Type==grid.LinearInterpFunction.EQUI
                % SSGrid is TensorGrid object in this case
                % reshape values for call to spapi
                dimvec=sf.SSGrid.Dimvec';
                ndim=sf.SSGrid.Ndim;
                vals=reshape(vals,[dimvec,nof]);
                vals=permute(vals,[ndim+1,1:ndim]);
                baseindex=repmat({':'},1,ndim+1);
                sf.InterpStruct=cell(nof,1);
                for f=1:nof
                    index=baseindex;
                    index{1}=f;
                    sf.InterpStruct{f}=griddedInterpolant(sf.SSGrid.Unigrids,squeeze(vals(index{:})),'linear','none');
                end
            else
                % SSGrid is ScatterGrid
                sf.InterpStruct=[];
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
            % check case
            if sf.Type==grid.LinearInterpFunction.EQUI
                vals=zeros(sf.Nof,np);
                for f=1:sf.Nof
                    vals(f,:)=sf.InterpStruct{f}(points_corr);
                end
            else
                % simplex search using saved triangulation
                tri=sf.SSGrid.Tessel;
                [Tin,bcc] = tsearchn(sf.SSGrid.Pointmat,tri,points_corr);
                K = ~isnan(Tin);
                vals = zeros(np,sf.Nof);
                for d = 1:ndim+1
                    % delaunay interp in each dimension
                    vals(K,:) = vals(K,:) + bsxfun(@times,bcc(K,d),sf.Vals(tri(Tin(K),d),:));
                end
                % for points outside of the convex hull of the standard
                % triangulation, use furthest-site triangulation instead
                % (option 'Qu' in Qhull)
                if any(~K)
                      % find closest point
                      [cpi1,d1]=dsearchn(sf.SSGrid.Pointmat,tri,points_corr(~K,:));
                      redind=setdiff(1:sf.SSGrid.Npt,cpi1);
                      % find second closest point
                      ptmat2=sf.SSGrid.Pointmat(redind,:);
                      fctvals2=sf.Vals(redind,:);
                      [cpi2,d2]=dsearchn(ptmat2,points_corr(~K,:));
                      % distance between closest and 2nd closest
                      b=(ptmat2(cpi2,:)-sf.SSGrid.Pointmat(cpi1,:)).^2;
                      b=sum(b,2);
                      % some simple trigonometry to get location on edge
                      % between to closest grid points that is overall closest point on
                      % convex hull
                      x=(d1.^2+b-d2.^2)./(2*b);
                      vals(~K,:)=sf.Vals(cpi1,:) + repmat(max(0,x),1,sf.Nof).*(fctvals2(cpi2,:)-sf.Vals(cpi1,:));
                end
                vals=vals';
            end
        end
        
        function plot2D(cf,dispfct,dispdim,val_other,ln)         
           
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
           figure;
           plot(gridx,vlist,'.-');
       end        
    
        
    end
    
    
end
