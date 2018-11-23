classdef dcimgReaderMatlab<handle
    properties
        metadata
        memmap
    end
% frame index 1 is the first frame. frameindex 0 does not exists in matlabs world
    methods
    function obj = dcimgReaderMatlab(fname)

        fId = fopen(fname, 'rb');
        % read from file
        fContent = fread(fId, 200, 'uint32'); % i.e. 4 byte
        fclose(fId);
        
        h = parseheader(fContent, fname);
        obj.metadata=h;
        obj.metadata.filename=fname;
        obj.memmap=memmapfile(fname, ...
            'Format', {'uint16', [double((h.bytes_per_frame / h.bytes_per_pixel)) h.num_frames], 'frame'},...
            'Offset', h.session0_data,'Repeat', 1);
    end

%     function frames_with_footer = getAllFrames(obj)
%     % functiont o read all frames
%         frames_with_footer = obj.memmap.Data.frame;
%     end

    function frames = getSpecificFrames(obj, frameIndices)
    % function to read only a subset of all frames
        h=obj.metadata;
        frames = zeros(h.num_columns,h.num_rows, length(frameIndices), 'uint16');
        for k = 1:numel(frameIndices)
            frames(:,:,k) = obj.getSingleReshapedFrame(obj.memmap.Data.frame(:,frameIndices(k)));
        end
    end

    function frame = getSingleReshapedFrame(obj,framesPlusFooter)
        h=obj.metadata;
        frames = framesPlusFooter(1:h.bytes_per_image / h.bytes_per_pixel);
        frame = (reshape(frames, h.num_columns, h.num_rows));
    end
    function close(obj)
        obj.memmap=[];
%         clear(obj.memmap)
    end
    end
end

    function h = parseheader(in, fname)

        % SESSION_HEADER_DTYPE_INT2P24

        h.format_version = in(3);

        if sum(ismember(h.format_version, [7 16777216])) <= 0
            warning("Format version is not 0x7 or  0x1000000 as expected");
        end

        h.num_sessions = in(9);
        if h.num_sessions ~= 1
            warning("It appears that there are %d sessions. We only support one.", h.num_sessions);
        end

        h.num_frames = in(10);
        h.session0_offset = in(11);
        h.filesize = in(13);
        h.filesize = in(17);
        h.mystery1 = in(22);

        % 
        if h.format_version == 7
            m = memmapfile(fname, ...
            'Format', { ...
                'uint32', 1, 'session_length';...
                'uint32', 1, 'pad1';...
                'uint32', 1, 'pseudo_offset';...
                'uint32', 5, 'pad2';...
                'uint32', 1, 'num_frames';...
                'uint32', 1, 'pixel_type';...
                'uint32', 1, 'mystery1';...
                'uint32', 1,'num_columns';...
                'uint32', 1,'bytes_per_row';...
                'uint32', 1,'num_rows';...
                'uint32', 1,'bytes_per_image';...
                'uint32', 2, 'pad3';...
                'uint32', 1,'offset_to_data';...
                'uint32', 1,'offset_to_footer' },...
                'Offset', h.session0_offset, ...
                'Repeat', 1);

        elseif  h.format_version == 16777216
            m = memmapfile(fname, ...
            'Format', { ...
                'uint32', 1,'session_length';...
                'uint32', 5, 'pad1';...
                'uint32', 1, 'pseudo_offset';...
                'uint32', 8, 'pad2';...
                'uint32', 1, 'num_frames';...
                'uint32', 1, 'pixel_type';...
                'uint32', 1, 'mystery1';...
                'uint32', 1, 'num_columns';...
                'uint32', 1, 'num_rows';...
                'uint32', 1, 'bytes_per_row';...
                'uint32', 1, 'bytes_per_image';...
                'uint32', 2, 'pad3';...
                'uint32', 1, 'offset_to_data';...
                'uint32', 4, 'pad4';...
                'uint32', 1, 'bytes_per_frame'},...
                'Offset', h.session0_offset, ...
                'Repeat', 1);
        end

        session_head = m.Data;

        h.pixel_type = double(session_head.pixel_type);
        h.num_columns = double(session_head.num_columns);
        h.bytes_per_row = double(session_head.bytes_per_row);
        h.bytes_per_pixel = double(h.bytes_per_row / h.num_columns);
        h.num_rows = double(session_head.num_rows);
        h.session0_data = double(h.session0_offset + session_head.offset_to_data);

        h.bytes_per_image = double(session_head.bytes_per_image);
        try
            h.bytes_per_frame = double(session_head.bytes_per_frame);
        catch
            % add ValueError catch here
%             warning("WARNING");
            h.bytes_per_frame = double(session_head.bytes_per_image);
        end
        clear m;

    end
    