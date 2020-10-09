classdef ScatterGrid < grid.StateSpaceGrid
    
    properties (SetAccess=protected)
        % inherited properties (abstract in superclass)
        Ndim
        Npt
        Pointmat
        Type
        % triangulation for scatter interp
        Tessel
    end
    
    properties (Constant)
        DelaunayOptions = {'Qt','Qbb','Qc','Qx'};
    end
    
    methods
        % constructor
        function ssg = ScatterGrid(arg1,arg2)
            if nargin < 1
                error('Not enough input arguments.');
            end
            if nargin==1
                % if initialized by specifying point matrix directly
                ssg.Pointmat=arg1;
                ssg.Type='scatter (direct)';
            else
                % if initialized by specifying tensor grid and index vector
                ptmat=arg1.Pointmat;
                imat=logical(arg2);
                if size(imat,1)~=size(ptmat,1)
                    error('index vector (arg2) must have same number of rows as point matrix of tensor grid (arg1)');
                end
                indvec=1:arg1.Npt;
                ssg.Pointmat=ptmat(indvec(imat),:);
                ssg.Type='scatter (tensor based) ';
            end
            ssg.Npt=size(ssg.Pointmat,1);
            ssg.Ndim=size(ssg.Pointmat,2);
            ssg.StateBounds=zeros(2,ssg.Ndim);
            for d=1:ssg.Ndim
                gridd=ssg.Pointmat(:,d);
                ssg.StateBounds(:,d)=[min(gridd);max(gridd)];
            end
            % precompute Delaunay Tessellation
            ssg.Tessel=delaunayn(ssg.Pointmat,grid.ScatterGrid.DelaunayOptions);
        end
           
        
        % add points to scatter grid
        function [ssg,newind,oldind] = addPoints(ssg,points)
            nptold=ssg.Npt;
            ptmat=[ssg.Pointmat; points];
%             [ptmat_un,~,ind_un]=unique(ptmat,'rows');
%             oldind=ind_un(1:nptold);
%             newind=ind_un(nptold+1:end);
            ptmat_un=ptmat;
            oldind=1:nptold;
            nptnew=size(points,1);
            newind=nptold+1:nptold+nptnew;
            ssg.Pointmat=ptmat_un;
            ssg.Npt=size(ssg.Pointmat,1);
            for d=1:ssg.Ndim
                gridd=ssg.Pointmat(:,d);
                ssg.StateBounds(:,d)=[min(gridd);max(gridd)];
            end          
            ssg.Tessel=delaunayn(ssg.Pointmat,grid.ScatterGrid.DelaunayOptions);
        end
        
        
    end
    
    
    
end