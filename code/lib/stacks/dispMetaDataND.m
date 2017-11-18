
fname = dir('*.nd2')
data = bfopen(fname.name);
metadata = data{1, 2};
metadataKeys = metadata.keySet().iterator();
for i=1:metadata.size()
  key = metadataKeys.nextElement();
  if strfind(key,'Global Z position for position, plane'), continue;end;
  if strfind(key,'Z position for position, plane'), continue;end;
  if strfind(key,'timestamp #'), continue;end;
  
  
  %if isempty(strfind(key,'Line:1;')), continue;end;
  %if isempty(strfind(key,'Equidistant')), continue;end;
  %if ~isempty(strfind(key, 'Z Stack Loop'))
      
      value = metadata.get(key);
      fprintf('%s = %s\n', key, value)
  %end
end
