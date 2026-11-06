function h=study_hash(value,isFile)
if nargin<2, isFile=false; end
if isFile
    f=fopen(value,'rb'); assert(f>=0); c=onCleanup(@()fclose(f)); %#ok<NASGU>
    bytes=fread(f,Inf,'*uint8');
else
    bytes=unicode2native(jsonencode(value),'UTF-8');
end
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(uint8(bytes));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
