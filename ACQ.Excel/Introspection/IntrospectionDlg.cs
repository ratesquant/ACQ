using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace ACQ.Excel.Introspection
{
    public partial class IntrospectionDlg : Form
    {
        ObjectIntrospector m_introspector;

        public IntrospectionDlg(object acq_object)
        {
            m_introspector = new ObjectIntrospector(acq_object);

            InitializeComponent();

            this.propertyGrid1.SelectedObject = m_introspector.SelectedObject;
            this.toolStripStatusLabel1.Text = m_introspector.Name;
            this.toolStripStatusLabel2.Text = String.Empty; //reserved

            if (m_introspector.IsDataObject)
            {
                this.dataGridView1.DataSource = m_introspector.Data;
                this.tabControl1.SelectedTab = this.tabPage2;
            }
            
        }
    }
}
