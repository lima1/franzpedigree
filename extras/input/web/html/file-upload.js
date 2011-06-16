Ext.onReady(function(){

    Ext.QuickTips.init();

    var msg = function(title, msg){
        Ext.Msg.show({
            title: title, 
            msg: msg,
            minWidth: 200,
            modal: true,
            icon: Ext.Msg.INFO,
            buttons: Ext.Msg.OK
        });
    };
        
    
    var tabs = new Ext.TabPanel({
        renderTo: 'csvresults',
        activeTab: 0,
        width: 590,
        height: 250,
        items: [ {
        id: 'ratab',   
        title: 'Main input file',
        xtype: 'textarea',
        },{
        id: 'pedtab',
        title: 'Pedigree file',
        xtype: 'textarea',
        }]
    });

    var fp = new Ext.FormPanel({
        renderTo: 'csvupper',
        fileUpload: true,
        width: 590,
        frame: true,
        title: 'FRANz CSV Importer 0.4.4',
        autoHeight: true,
        bodyStyle: 'padding: 10px 10px 0 10px;',
        labelWidth: 190,
        defaults: {
            anchor: '95%',
            allowBlank: false,
            msgTarget: 'side'
        },
        items: [
        {
            xtype: 'checkbox',
            id: 'has_header',
            checked: true,
            fieldLabel: 'Has Header Row',
        },
        {
            xtype: 'checkbox',
            id: 'read_locus_names',
            checked: true,
            fieldLabel: 'Read Locus Names',
        },
        {
            xtype: 'textfield',
            id: 'sep_char',
            fieldLabel: 'Separator',
            value: ',',
            maxLength: 1,
        },
        {
            xtype: 'textfield',
            id: 'dataset_title',
            fieldLabel: 'Data-set Title',
            value: 'FRANz',
        },
        {
            xtype: 'numberfield',
            id: 'id_col',
            fieldLabel: 'Genotype IDs in column',
            value: 1,
        },
        {
            xtype: 'numberfield',
            id: 'sex_col',
            fieldLabel: 'Sex in column',
            emptyText: 'optional',
            allowBlank: true,
        },
        {
            xtype: 'numberfield',
            id: 'birth_col',
            fieldLabel: 'Year of birth in column',
            emptyText: 'optional',
            allowBlank: true,
        },
        {
            xtype: 'numberfield',
            id: 'death_col',
            fieldLabel: 'Year of death in column',
            emptyText: 'optional',
            allowBlank: true,
        },
        {
            xtype: 'numberfield',
            id: 'mother_col',
            fieldLabel: 'ID of mother in column',
            emptyText: 'optional',
            allowBlank: true,
        },
        {
            xtype: 'numberfield',
            id: 'location_col',
            fieldLabel: 'ID of sampling location in column',
            emptyText: 'optional',
            allowBlank: true,
        },
        {
            xtype: 'numberfield',
            id: 'data_col',
            fieldLabel: 'Alleles start at column',
            value: 2,
        },
        {
            xtype: 'fileuploadfield',
            id: 'form-file',
            emptyText: 'Select a CSV file',
            fieldLabel: 'CSV',
            name: 'csv-path',
            buttonCfg: {
                text: 'Browse',
            }
        }],
        buttons: [{
            text: 'Save',
            handler: function(){
                if(fp.getForm().isValid()){
                    fp.getForm().submit({
    url: 'http://www.bioinf.uni-leipzig.de/cgi-bin/markus/file-upload.cgi',
                        waitMsg: 'Uploading your data...',
                        success: function(fp, o){
                            msg('Success', 'Processed file. If everything looks good, paste the content below in a new text file! If not, change the settings and hit the Save-Button again.' );
        Ext.getCmp('ratab').setRawValue(o.result.data);
        Ext.getCmp('pedtab').setRawValue(o.result.ped);
                            //ra.show();
                        },
                        failure: function(fp, o){
                            msg('Failure', o.result.msg );
                        }
                    });
                }
            }
        },{
            text: 'Reset',
            handler: function(){
                fp.getForm().reset();
            }
        }]
    });


});
