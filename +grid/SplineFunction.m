classdef SplineFunction < grid.ApproxFunction
    
    properties (SetAccess=protected)
        % inherited properties (abstract in superclass)
        SSGrid
        Nof
        Vals      
        % spline specific properties
        Degvec
        SplineStruct
        Extrap
    end
    
    methods
        % constructor
        function sf=SplineFunction(ssgrid,vals,degvec,extrap)
            sf.SSGrid=ssgrid;
            sf.Degvec=degvec;
            sf.Extrap=extrap;
            sf=fitTo(sf,vals);
        end
        
        function sf=set.SSGrid(sf,ssg)
            % check that grid object at initialization is TensorGrid
            if ~isa(ssg,'grid.TensorGrid')
                error('StateSpaceGrid must be a TensorGrid');
            end
            sf.SSGrid=ssg;
        end
        
        function sf=set.Degvec(sf,degvec)
            if size(degvec,1)>1
                degvec=permute(degvec,[2,1]);
            end
            if size(degvec,1)>1
                error('Degvec must be a vector of length Ndim');
            end
            sf.Degvec=degvec;
        end
        
        function sf=set.Extrap(sf,extrap)
            if size(extrap,1)>2
                extrap=permute(extrap,[2,1]);
            end
            if size(extrap,1)>2
                error('Extrap must have dimension 2 x Ndim');
            end
            sf.Extrap=extrap;
        end
        
        % fit to new values
        function sf=fitTo(sf,vals)
            [npt,nof]=size(vals);
            if npt~=sf.SSGrid.Npt
                error('Value matrix must have dimensions (Npt x Nof)');
            end
            sf.Nof=nof;
            sf.Vals=vals;
            % get some basics from grid object
            dimvec=sf.SSGrid.Dimvec';
            ndim=sf.SSGrid.Ndim;
            % reshape values for call to spapi
            vals=reshape(vals,[dimvec,nof]);
            vals=permute(vals,[ndim+1,1:ndim]);
            % call to spapi from matlab spline toolbox
            if ndim==1
                % workaround for bug in fnxtr
                spl=spapi(sf.Degvec,sf.SSGrid.Unigrids{1},vals);
            else
                spl=spapi(num2cell(sf.Degvec),sf.SSGrid.Unigrids,vals);                
            end
            % make pp form for ppmval
            sf.SplineStruct=fn2fm(spl,'pp');
        end
        
        % evaluation
        function vals=evaluateAt(sf,points)
            [np,ndim]=size(points);
            if ndim~=sf.SSGrid.Ndim
                error('Point matrix must have dimensions (#points x Ndim)');
            end
            % check extrapolation rules
            fbnds_down=ones(np,1)*sf.Extrap(1,:);
            fbnds_up=ones(np,1)*sf.Extrap(2,:);
            SBlow=ones(np,1)*sf.SSGrid.StateBounds(1,:);
            SBhi=ones(np,1)*sf.SSGrid.StateBounds(2,:);
            points_corr=points;
            upvio=logical(fbnds_up.*(points>SBhi));
            points_corr(upvio)=SBhi(upvio);
            downvio=logical(fbnds_down.*(points<SBlow));
            points_corr(downvio)=SBlow(downvio);
            % call to ppmval (NOTE: currently only works up to 4 dim)
            try
                % MEX function
                vals=ppmval(points_corr',sf.SplineStruct);
            catch
                % Matlab function (for non-Windows systems)
                vals=fnval(sf.SplineStruct,points_corr');
            end
        end
        
        
    end
    
    
end
