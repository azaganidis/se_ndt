%plots and processes the covariance log file

[n_scans n_samples n_dim] = size(icp);

for i=1:n_scans
    figure;
    for j1=1:n_dim
        for j2=j1+1:n_dim
            subplot(n_dim-1, n_dim-1, (j1-1)*(n_dim-1)+j2-1);
            hold on;
            %TODO: go through and remove outliers... they only screw up the
            %data
            plot(icp(i,:,j1),icp(i,:,j2),'rx', ...
                pnt2ndt(i,:,j1),pnt2ndt(i,:,j2),'ks', ...
                ndt2ndt(i,:,j1),ndt2ndt(i,:,j2),'bo');
            Cicp = cov([icp(i,:,j1);icp(i,:,j2)]');
            micp = [mean(icp(i,:,j1));mean(icp(i,:,j2))];
            Cpnt2ndt = cov([pnt2ndt(i,:,j1);pnt2ndt(i,:,j2)]');
            mpnt2ndt = [mean(pnt2ndt(i,:,j1));mean(pnt2ndt(i,:,j2))];
            Cndt2ndt = cov([ndt2ndt(i,:,j1);ndt2ndt(i,:,j2)]');
            mndt2ndt = [mean(ndt2ndt(i,:,j1));mean(ndt2ndt(i,:,j2))];
            
            error_ellipse(Cicp,micp,'style','r');
            error_ellipse(Cpnt2ndt,mpnt2ndt,'style','k');
            error_ellipse(Cndt2ndt,mndt2ndt,'style','b');
            axis square;
            
            %title('error plot in dimensions '+j1+'-'+j2);
        end
    end
    
end
