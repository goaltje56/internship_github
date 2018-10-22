function [El,Sp,Re]=ReadTrotDat(filename, varargin)
% Function ReadThermDat reads properties from ALEX chemistry file
%  Output: El
%              El.Name
%              El.Mass
%          Sp
%              Sp.Name
%              Sp.Mass
%              Sp.pol  [2x8]     JANAF polynimials for Thermodynamic data
%              Sp.visc [var]     Polnomial coefficients for viscosity
%              Sp.viscord        Order of viscosity polynomial
%              Sp.cond [var]     Polnomial coefficients for conductivity
%              Sp.condord        Order of conductivity polynomial
%              Sp.diff [NspxNsp] Polynomial coefficients for binary 
%                                diffusion coefficients
%          Both structs!
clear Sp El Re

if nargin > 2
    error('Only one optional argument possible')
elseif nargin == 2
    show = varargin{1};
else
    show = true;
end

fid=fopen(filename);
SFP=ftell(fid);
Ready=0;
Keys(1).Name='ELEM';
Keys(2).Name='SPEC';
Keys(3).Name='THERMO';
Keys(4).Name='VISCOSITY';
Keys(5).Name='CONDUCTIVITY';
Keys(6).Name='DIFFUSIVITIES';
Keys(7).Name='REAC';
Found=0;iKey=1;NumKeys=length(Keys);
while (iKey<=NumKeys)
    fseek(fid,SFP,'bof');
    CurKey=Keys(iKey).Name;
    KeyLength=length(CurKey);
    while ((Found==0)&(feof(fid)==0))
        CFP=ftell(fid);
        CurLine=fgetl(fid);LineLength=length(CurLine);
        Lend=min([KeyLength+5 LineLength]);
        Matches=findstr(CurKey,CurLine(1:Lend));
        if (length(Matches)>0)&(LineLength>KeyLength)
            Found=1;
            iKey=iKey+1;
        end;
        clear Matches;
    end;
    if (Found==1)
        fseek(fid,CFP,'bof');
        switch CurKey
            case 'ELEM'
                Check=fscanf(fid,'%s',[1,1]);
                Ne=fscanf(fid,'%i',[1,1]);
                if show
                    disp(['Processing [' Check ']: ' num2str(Ne) ' entries'])
                end
                for i=1:Ne
                    El(i).Name=fscanf(fid,'%s',[1,1]);
                    El(i).Mass=fscanf(fid,'%f',[1,1]);
                end;
            case 'SPEC'
                Check=fscanf(fid,'%s',[1,1]);
                Nsp=fscanf(fid,'%i',[1,1]);
                if show
                    disp(['Processing [' Check ']: ' num2str(Nsp) ' entries'])
                end
                for i=1:Nsp
                    Sp(i).Name=fscanf(fid,'%s',[1,1]);
                end;
            case 'REAC'
                Check=fscanf(fid,'%s',[1,1]);
                Nreac=fscanf(fid,'%i',[1,1]);
                if show
                    disp(['Processing [' Check ']: ' num2str(Nreac) ' entries'])
                end
                Diaz=fscanf(fid,'%s',[1,1]);    % #
                for k=1:Nreac
                    ireac=fscanf(fid,'%i',[1,1]); % number of reaction
                    ifile=fscanf(fid,'%i',[1,1]); % number of file
                    iline=fscanf(fid,'%i',[1,1]); % number of line in file
                    rev=fscanf(fid,'%s',[1,1]);   % Reversible
                    if strcmp(rev,'REV')
                        Re(k).rev=true;
                    elseif strcmp(rev,'IRREV')
                        Re(k).rev=false;
                    elseif show
                        disp('Unknown value');
                    end
                    pdep=fscanf(fid,'%s',[1,1]);  % Pressure dependent
                    if strcmp(pdep,'PDEP')
                        Re(k).pdep=true;
                    elseif strcmp(pdep,'INDEP')
                        Re(k).pdep=false;
                    elseif show
                        disp('Unknown value');
                    end
                    
                    third=fscanf(fid,'%s',[1,1]); % Third body
                    
                    ford=fscanf(fid,'%i',[1,1]);  % Number of species lhs
                    for i=1:ford
                        nu=fscanf(fid,'%i',[1,1]);  % Stoiciometric coefficient
                        pow=fscanf(fid,'%i',[1,1]); % Power
                        Name=fscanf(fid,'%s',[1,1]); % Species name
                        Re(k).fspec(i)=find(strcmp({Sp.Name,'M'},Name));
                        Re(k).fstoi(i)=nu;
                        Re(k).fpow(i)=pow;
                    end
                    
                    rord=fscanf(fid,'%i',[1,1]);  % Number of species rhs
                    for i=1:rord
                        nu=fscanf(fid,'%i',[1,1]);  % Stoiciometric coefficient
                        pow=fscanf(fid,'%i',[1,1]); % Power
                        Name=fscanf(fid,'%s',[1,1]); % Species name
                        Re(k).rspec(i)=find(strcmp({Sp.Name,'M'},Name));
                        Re(k).rstoi(i)=nu;
                        Re(k).rpow(i)=pow;
                    end
                    
                    % Arrhenius coefficients
                    Re(k).preexp=fscanf(fid,'%f',[1,1]);
                    Re(k).Tpower=fscanf(fid,'%f',[1,1]);
                    Re(k).Tactivation=fscanf(fid,'%f',[1,1]);
                    
                    Check=fscanf(fid,'%s',[1,1]);
                    if strcmp(Check,Diaz)||strcmp(Check,'END')
                        auxinfo=false; 
                    else
                        auxinfo=true;
                    end
                    
                    Re(k).revgiven=false;
                    Re(k).revdata=[];
                    Re(k).enhanced=false;
                    Re(k).enhspec=[];
                    Re(k).enhval=[];
                    Re(k).low=false;
                    Re(k).lowdata=[];
                    Re(k).troe=false;
                    Re(k).troedata=[];
                    Re(k).sri=false;
                    Re(k).sridata=[];
                    Re(k).dup=false;
                    
                    while auxinfo
                        if strcmp(Check,'REV')
                            Re(k).revgiven=true;
                            Re(k).revdata=fscanf(fid,'%f',[1,3]);
                        elseif strcmp(Check,'ENHANCED')
                            Re(k).enhanced=true;
                            Nenh=fscanf(fid,'%i',[1,1]);
                            for i=1:Nenh
                                Name=fscanf(fid,'%s',[1,1]); % species name
                                val=fscanf(fid,'%f',[1,1]); % value
                                Re(k).enhspec(i)=find(strcmp({Sp.Name},Name));
                                Re(k).enhval(i)=val;
                            end
                        elseif strcmp(Check,'LOW')
                            Re(k).low=true;
                            Re(k).lowdata=fscanf(fid,'%f',[1,3]);
                        elseif strcmp(Check,'TROE')
                            Re(k).troe=true;
                            Ntroe=fscanf(fid,'%i',[1,1]);
                            Re(k).troedata=fscanf(fid,'%f',[1,Ntroe]);
                        elseif strcmp(Check,'SRI')
                            Re(k).sri=true;
                            Re(k).sridata=fscanf(fid,'%f',[1,6]);
                        elseif strcmp(Check,'DUP')
                            Re(k).dup=true;
                        elseif show
                            disp('Unknown auxiliary information');
                        end

                        Check=fscanf(fid,'%s',[1,1]);
                        if strcmp(Check,Diaz)||strcmp(Check,'END')
                            auxinfo=false;
                        else
                            auxinfo=true;
                        end
                    end
                end;
            case 'THERMO'
                Check=fscanf(fid,'%s',[1,1]);
                Nsp=fscanf(fid,'%i',[1,1]);
                if show
                    disp(['Processing [' Check ']: ' num2str(Nsp) ' entries'])
                end
                for ii=1:Nsp
                    SkipLine=fgetl(fid);
                    Name=fscanf(fid,'%s\n',[1,1]);
                    Phase=fscanf(fid,'%s\n',[1,1]);
                    %                     SkipLine=fgetl(fid);
                    CurNel=fscanf(fid,'%i\n',[1,1]);
                    %                     SkipLine=fgetl(fid);
                    comp=zeros(size([El.Mass]));
                    for i=1:CurNel
                        cdum=fscanf(fid','%s',[1,1]);
                        cnum=fscanf(fid,'%i');
                        for j=1:Ne,
                            if (El(j).Name==cdum)
                                comp(j)=cnum;
                            end;
                        end;
                    end;
                    dummy=fscanf(fid,'%s',[1,1]);
                    Sp(ii).Mass=fscanf(fid,'%f\n',[1,1]);
                    SkipLine=fgetl(fid);
                    T1range=fscanf(fid,'%f',[1,2]);
                    Cp1=fscanf(fid,'%f %f',[1,8]);
                    T2range=fscanf(fid,'%f',[1,2]);
                    Cp2=fscanf(fid,'%f',[1,8]);
                    Sp(ii).Ts=T1range(2);
                    Sp(ii).pol=[Cp1;Cp2];
                    Sp(ii).elcomp=comp;
                    Sp(ii).phase=Phase;
                    Sp(ii).Comment=fgetl(fid);
                end;
            case 'VISCOSITY'
                Check=fscanf(fid,'%s',[1,1]);
                Nsp=fscanf(fid,'%i',[1,1]);
                if show
                    disp(['Processing [' Check ']: ' num2str(Nsp) ' entries'])
                end
                for ii=1:Nsp
                    Name=fscanf(fid,'%s',[1,1]);
                    Ord=fscanf(fid,'%i',[1,1]);
                    Sp(ii).viscord=Ord;
                    Sp(ii).visc=fscanf(fid,'%f',[1,Ord]);
                end;
            case 'CONDUCTIVITY'
                Check=fscanf(fid,'%s',[1,1]);
                Nsp=fscanf(fid,'%i',[1,1]);
                if show
                    disp(['Processing [' Check ']: ' num2str(Nsp) ' entries'])
                end
                for ii=1:Nsp
                    Name=fscanf(fid,'%s',[1,1]);
                    Ord=fscanf(fid,'%i',[1,1]);
                    Sp(ii).condord=Ord;
                    Sp(ii).cond=fscanf(fid,'%f',[1,Ord]);
                end;
            case 'DIFFUSIVITIES'
                Check=fscanf(fid,'%s',[1,1]);
                Nsp=fscanf(fid,'%i',[1,1]);
                if show
                    disp(['Processing [' Check ']: ' num2str(Nsp) ' entries'])
                end
                for ii=1:Nsp
                    for jj=1:ii
                        Name1=fscanf(fid,'%s',[1,1]);
                        Name2=fscanf(fid,'%s',[1,1]);
                        Ord=fscanf(fid,'%i',[1,1]);
                        Sp(ii).difford(jj)=Ord;
                        Sp(ii).diff(jj,:)=zeros(1,10);
                        Sp(ii).diff(jj,1:Ord)=fscanf(fid,'%f',[1,Ord]);
                        %fscanf(fid,'%f',[1,Ord]);
                    end;
                end;
%               Mirror the matrix D_ij=D_ji                
                for ii=1:Nsp
                    for jj=ii+1:Nsp
                        Sp(ii).difford(jj)=Sp(jj).difford(ii);
                        Sp(ii).diff(jj,:)=Sp(jj).diff(ii,:);
                    end;
                end;
            otherwise
                if show
                    disp(['Yet undefined case: [' CurKey ']'])
                end
        end;
        Found=0;
    else
        if show
            disp([CurKey ' NOT found'])
            disp(['Apparently invalid!'])
        end
        iKey=iKey+1;
    end;
end;
% clear Check Ne Nsp CurKey CurLine Name Ord Cp1 Cp2 T1range T2range Ready Found
% clear CFP SFP NumKeys
fclose(fid);
