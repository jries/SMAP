    function    str  = current_class( )
        dbk = dbstack( 1, '-completenames' );
        if isempty( dbk )
            str = 'base';
        else
            f=dbk(1).file;
            [~,str]=fileparts(f);
%             str = dbk(1).name;
        end
%         cac = regexp( str, '\.', 'split' );
%         switch  numel( cac )
%             case 1
%             case 2
%                 cls = meta.class.fromName( cac{1} );
%                 str = cac{2};
%             case 3
%             otherwise
%         end
    end