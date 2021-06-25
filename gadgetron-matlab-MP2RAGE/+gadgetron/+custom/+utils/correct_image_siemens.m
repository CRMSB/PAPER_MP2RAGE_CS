function image = correct_image_siemens(image,acquisition, connection)
headerScan=structfun(@(arr) arr(:, end)', acquisition.data.header, 'UniformOutput', false);


% prepare parameters
fov = [connection.header.encoding.reconSpace.fieldOfView_mm.x...
    connection.header.encoding.reconSpace.fieldOfView_mm.y...
    connection.header.encoding.reconSpace.fieldOfView_mm.z
    ];
matrixSize = [connection.header.encoding.reconSpace.matrixSize.x...
    connection.header.encoding.reconSpace.matrixSize.y...
    connection.header.encoding.reconSpace.matrixSize.z
    ];

res = fov./matrixSize;

read_dir = headerScan.read_dir;
phase_dir = headerScan.phase_dir;
slice_dir = headerScan.slice_dir;

vec = [-res(1); -res(2); -res(3)/2];
Mrot = [read_dir', phase_dir', slice_dir'];

% replace header fields
image.data = single(image.data);
image.header.position = headerScan.position' + Mrot * vec;
image.header.field_of_view = fov;
image.header.repetition = headerScan.repetition;
image.header.image_index = headerScan.repetition + 1;
end


