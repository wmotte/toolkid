% Name : procmshotepibl
%
% Purpose : To reconstruct a series of multishot multislice EPI images
% obtained in a series of repetititions which need to be averaged before
% Fourier transformation. 
%
% Calling sequence : im=procmshotepibl(refname,imname,exportname,zffactor1);
%
% Input : refname : Name of file containing the reference images
%         imname : Name of file containing the EPI series
%         exportname : Name of .bfloat file in which the images will be
%               exported.
%         zffactor1 : Zero fill multiplication factor in read-out direction
%
% Output : The images in .bfloat format arranged per slice
%
% Essential subroutines : read_epifid
%                         read_procpar
%                         epiphasemap
%                         phasewithmap
%                         saveblock
%                         remdcoffset
%
% Modification history : Author : Annette van der Toorn 6-1-2005

function im=procmshotepilc(refname,imname,expname,zffactor1,varargin);

startdir=cd;

if nargin > 4
   type=varargin{1};
else
   type='bfloat';
end

[np,nv,ntraces,nblocks]=read_epifidbl(refname);

zf1=pow2(ceil(log2(np-2)))*zffactor1;
zf2=zf1;
ncoils=nblocks;

pss=read_procpar(refname,'pss');
[pssord,pssindex]=sort(pss);

ns=length(pss);

phasemap=zeros(zf1,nv,ntraces,nblocks);

for blockcount=1:nblocks;
    tempfilename=strcat(refname,'.fid/mat',num2str(blockcount));
    load(tempfilename,'fid');
    for tracecount=1:ns;
        phasemap(:,:,tracecount,blockcount)=epiphasemap(squeeze(fid(:,:,pssindex(tracecount))),zf1);
    end;
end;

tempfilename=strcat(refname,'.fid/phasemap');
save(tempfilename,'phasemap');

[np,nv,ntraces,nblocks]=read_epifidbl(imname);

zf1=pow2(ceil(log2(np-2)))*zffactor1;
zf2=zf1;

pss=read_procpar(refname,'pss');
ns=length(pss);

pss=read_procpar(imname,'pss');
ns=length(pss);
seqcon=read_procpar(imname,'seqcon');
korder=read_procpar(imname,'ky_order');

[pssord,pssindex]=sort(pss);

nv2=read_procpar(imname,'nv2');

switch seqcon;
    case 'ccncn'; % look-locker compressed
        if nv2 >1;
            %nsteps=ntraces/(ns*nv2);
            
            trimage = read_procpar( imname, 'imagetime' );
            ti = read_procpar( imname, 'ti' );
            
            for i = 1:ns * nv2;
                titable( i ) = ti + ( i - 1 ) * trimage;
            end;

            titable = reshape( titable, ns, nv2 );
            
            for i = 1:ns;
                titemp = squeeze( titable( pssindex( i ), : ) );
            end;
             
                %newname=strcat(expname,'_ti','_',num2str(i-1));
                %fileid=fopen(newname,'w');
                %for j=1:size(titemp(:));
                %    fprintf(fileid,'%14.7f\n',titemp(j));
                %end;
                %fclose(fileid);
            end;


        else
            nsteps=ntraces/ns;
            nv2=1;
        end;
    case 'ccnsn'; % fMRI mode
        nsteps=ntraces/ns;
    case 'csncn'; % look-locker uncompressed slices
        if nv2>1;
            nsteps=ntraces/nv2;
            trimage=read_procpar(imname,'imagetime');
            ti=read_procpar(imname,'ti');
            for i=1:nv2;titemp(i)=ti+(i-1)*trimage;end;
            for i=1:ns;
                newname=strcat(expname,'_ti','_',num2str(i-1));
                fileid=fopen(newname,'w');
                for j=1:size(titemp(:));
                    fprintf(fileid,'%14.7f\n',titemp(j));
                end;
                fclose(fileid);
            end;
        else
            nsteps=ntraces;
            nv2=1;
        end;            
    otherwise;
        if nv2 >1;
            nsteps=ntraces/(ns*nv2);
            trimage=read_procpar(imname,'imagetime');
            ti=read_procpar(imname,'ti');
            for i=1:ns*nv2;titable(i)=ti+(i-1)*trimage;end;
            titable=reshape(titable,ns,nv2);
            for i=1:ns;
                titemp=squeeze(titable(pssindex(i),:));
                newname=strcat(expname,'_ti','_',num2str(i-1));
                fileid=fopen(newname,'w');
                for j=1:size(titemp(:));
                    fprintf(fileid,'%14.7f\n',titemp(j));
                end;
                fclose(fileid);
            end;
        else
            nsteps=ntraces/ns;
            nv2=1;
        end;        
end;

switch korder
    case 'c';
        t1=ones(nsteps,1)';
        t1(1:nsteps/2)=-1;
        t2=[-nsteps/2:1:nsteps/2-1];
        ktraj=zeros(nsteps,nv);
        for i=1:nsteps; 
            for j=1:nv;
                ktraj(i,j)=t2(i)+t1(i)*(j-1)*nsteps/2;
            end;
        end;
        ktraj=reshape(ktraj',nv*nsteps,1)
        [peord,peindex]=sort(ktraj);
    case 'l'
        ktraj=zeros(nsteps,nv);
        for i=1:nsteps;
            for j=1:nv;
                ktraj(i,j)=-nsteps/2*nv+(i-1)+nsteps*(j-1);
            end;
        end;
        ktraj=reshape(ktraj',nv*nsteps,1);
        [peord,peindex]=sort(ktraj);
    otherwise
        display('This ordering is not implemented');
end;
        

switch seqcon(2);
    case 'c';
        switch seqcon(4);
            case 'c';
                for blockcount=1:nblocks;
                    blockcount
                    coilcount=rem(blockcount-1,ncoils)+1;
                    tempfilename=strcat(imname,'.fid/mat',num2str(blockcount));
                    load(tempfilename,'fid');
                    fid=reshape(fid,np,nv,ns,nv2,nsteps);
                    for imagecount=1:nv2;
                        for slicecount=1:ns;
                            slicecount
                            mshotim=zeros(zf1,nv,nsteps);
                            fidtrace=squeeze(fid(:,:,pssindex(slicecount),imagecount,:));
                            for stepcount=1:nsteps;
                                tempim=fft1d(squeeze(fidtrace(:,:,stepcount)),zf1);
                                tempim=phasewithmap(tempim,phasemap(:,:,slicecount,coilcount));
                                mshotim(:,:,stepcount)=tempim;
                            end;
                            mshotim=reshape(mshotim,zf1,nv*nsteps);
                            tempim=mshotim;
                            for i=1:nv*nsteps
                                mshotim(:,i)=tempim(:,peindex(i));
                            end;
                            mshotim=rot90(mshotim,1);
                            mshotim=rot90(fft1d(mshotim,zf2),-1);
                            mshotim=flipdim(mshotim,2);
                            im(:,:,slicecount,imagecount,blockcount)=mshotim;
                        end;
                    end;
                end;
            case 's';
                for blockcount=1:nblocks;
                    blockcount
                    coilcount=rem(blockcount-1,ncoils)+1;
                    tempfilename=strcat(imname,'.fid/mat',num2str(blockcount));
                    load(tempfilename,'fid');
                    fid=reshape(fid,np,nv,ns,nsteps);
                    for slicecount=1:ns;
                        slicecount
                        mshotim=zeros(zf1,nv,nsteps);
                        fidtrace=squeeze(fid(:,:,pssindex(slicecount),:));
                        for stepcount=1:nsteps;
                            tempim=fft1d(squeeze(fidtrace(:,:,stepcount)),zf1);
                            tempim=phasewithmap(tempim,phasemap(:,:,slicecount,coilcount));
                            mshotim(:,:,stepcount)=tempim;
                        end;
                        mshotim=reshape(mshotim,zf1,nv*nsteps);
                        tempim=mshotim;
                        for i=1:nv*nsteps
                            mshotim(:,i)=tempim(:,peindex(i));
                        end;
                        mshotim=rot90(mshotim,1);
                        mshotim=rot90(fft1d(mshotim,zf2),-1);
                        mshotim=flipdim(mshotim,2);
                        im(:,:,slicecount,blockcount)=mshotim;
                    end;
                end;
            case 'n';
                im=zeros(zf1,zf2,ns,nblocks);
                for blockcount=1:nblocks;
                    blockcount
                    coilcount=rem(blockcount-1,ncoils)+1;
                    tempfilename=strcat(imname,'.fid/mat',num2str(blockcount));
                    load(tempfilename,'fid');
                    fid=reshape(fid,np,nv,ns,nsteps);
                    for slicecount=1:ns;
                        slicecount
                        mshotim=zeros(zf1,nv,nsteps);
                        fidtrace=squeeze(fid(:,:,pssindex(slicecount),:));
                        for stepcount=1:nsteps;
                            tempim=fft1d(squeeze(fidtrace(:,:,stepcount)),zf1);
                            tempim=phasewithmap(tempim,phasemap(:,:,slicecount,coilcount));
                            mshotim(:,:,stepcount)=tempim;
                        end;
                        mshotim=reshape(mshotim,zf1,nv*nsteps);
                        tempim=mshotim;
                        for i=1:nv*nsteps
                            mshotim(:,i)=tempim(:,peindex(i));
                        end;
                        mshotim=rot90(mshotim,1);
                        mshotim=rot90(fft1d(mshotim,zf2),-1);
                        mshotim=flipdim(mshotim,2);
                        im(:,:,slicecount,blockcount)=mshotim;
                    end;
                end;
        end;
    case 's';
        im=zeros(zf1,zf2,nv2,nblocks);
        for blockcount=1:nblocks; % Is also equal to nslices unless ncoils>1
            blockcount;
            coilcount=rem(blockcount-1,ncoils)+1;
            if ncoils > 1;
                slicecount=pssindex(floor(blockcount/ncoils)+1)
            else
                slicecount=pssindex(blockcount)
            end;
            tempfilename=strcat(imname,'.fid/mat',num2str(blockcount));
            load(tempfilename,'fid');
            fid=reshape(fid,np,nv,nv2,nsteps);
            for imagecount=1:nv2;
                mshotim=zeros(zf1,nv,nsteps);
                fidtrace=squeeze(fid(:,:,imagecount,:));
                for stepcount=1:nsteps;
                    tempim=fft1d(squeeze(fidtrace(:,:,stepcount)),zf1);
                    tempim=phasewithmap(tempim,phasemap(:,:,slicecount,coilcount));
                    mshotim(:,:,stepcount)=tempim;
                end;
                mshotim=reshape(mshotim,zf1,nv*nsteps);
                tempim=mshotim;
                for i=1:nv*nsteps
                    mshotim(:,i)=tempim(:,peindex(i));
                end;
                mshotim=rot90(mshotim,1);
                mshotim=rot90(fft1d(mshotim,zf2),-1);
                mshotim=flipdim(mshotim,2);
                im(:,:,imagecount,blockcount)=mshotim;
            end;
        end;
        im=reshape(im,zf1,zf2,nv2,nblocks/ns,ns);
        for slicecount=1:ns;
            imtemp(:,:,:,:,slicecount)=im(:,:,:,:,pssindex(slicecount));
        end;
        im=imtemp;
        im=permute(im,[1,2,5,3,4]);
    otherwise;
        display('Slices are always arranged somehow');
end;

im=squeeze(im);
szim=size(im);
if length(szim) < 5;
    if length(szim) < 4;
        if length(szim) < 3;
            szim(3)=1;
        end;
        szim(4)=1;
    end;
    szim(5)=1;
end;

im=reshape(im,szim(1),szim(2),szim(3),szim(4),ncoils,szim(5)/ncoils);
im=permute(im,[1,2,3,4,6,5]);
tempim=zeros(szim(1),szim(2),szim(3),szim(4),szim(5)/ncoils);
if coilcount >1;
    for coilcount=1:ncoils;
        tempim=tempim+abs(im(:,:,:,:,:,coilcount));
    end;
    im=tempim;
end;

switch type;
    case 'bfloat';
        exportbfloat(abs(im),expname);
    case 'sdt';
        save_sdt(abs(im),expname);
    case 'img';
        szim=size(im);
        if length(szim) > 3;
            for timecount=1:szim(4);
                strcount=num2str(1000+timecount);
                newname=strcat(expname,'_',strcount(2:4));
                save_img(squeeze(abs(im(:,:,:,timecount))),imname,newname);
            end;
        else;
            save_img(squeeze(abs(im)),imname,expname);
        end;
    case 'nii';
        save_nifti(abs(im),imname,expname);
    otherwise;
        disp('data type unknown');
end;


cd(startdir);
newname=strcat(imname,'.fid');
cd(newname);
delete('*.mat');
cd(startdir);
newname=strcat(refname,'.fid');
cd(newname);
delete('*.mat');
cd(startdir);


