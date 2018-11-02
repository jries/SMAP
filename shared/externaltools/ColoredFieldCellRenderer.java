// ColoredFieldCellRenderer - Modified TableCellRenderer for colored cells

// Programmed by Yair M. Altman: altmany(at)gmail.com

import java.awt.*;
import java.util.*;
import javax.swing.*;
import javax.swing.table.*;

public class ColoredFieldCellRenderer extends DefaultTableCellRenderer implements TableCellRenderer
{
	private Color _bgcolor = getBackground();
	private Color _fgcolor = getForeground();
	private boolean _debug = false;
	private boolean _disabled = false;
	private boolean _smartAlign = true;
	private Hashtable _cellBgColorHashtable = new Hashtable();
	private Hashtable _cellFgColorHashtable = new Hashtable();
	private Hashtable _cellTooltipHashtable = new Hashtable();

	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
	{
		JComponent cell = (JComponent) super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
		if (_debug) System.out.println(row + "," + column + " => value: " + value);

		// Use the chosen default bgcolor (possibly overridden later)
        String valueStr = "" + value;
		cell.setBackground(_bgcolor);
		cell.setForeground(_fgcolor);

        /* Try to fix sorted-row reordering
        try
        {
            TableSorter model = (TableSorter) table.getModel();
            row = model.modelIndex(row);
        }
        catch (Exception ex)
        {
            // never mind...
        }
        */

		if (_disabled)
		{
			cell.setEnabled(false);
			cell.setBackground(_bgcolor);
		}

		// Align the cell contents
		try
		{
			// If SmartAlign was requested
			if (_smartAlign)
			{
				// If this is a string, align left - otherwise right
				java.lang.String text = (java.lang.String) value;
				if  ((text.charAt(0) == '-') || (text.charAt(0) == '.') || 
					((text.charAt(0) >= '0') && (text.charAt(0) <= '9')))
				{
					((JLabel)cell).setHorizontalAlignment(SwingConstants.RIGHT);
				}
				else
				{
					((JLabel)cell).setHorizontalAlignment(SwingConstants.LEFT);
				}
			} else {
				// Use the CellRenderer's default HorizontalAlignment property
				((JLabel)cell).setHorizontalAlignment(getHorizontalAlignment());
			}
        }
		catch (Exception ex)
		{
			((JLabel)cell).setHorizontalAlignment(SwingConstants.RIGHT);
		}

		// If this cell should have a specific color, then use it
		Vector rowColVector = getRowColVector(row, column);
		Color cellBgColor = (Color) _cellBgColorHashtable.get(rowColVector);
		if (cellBgColor == null)
			cellBgColor = (Color) _cellBgColorHashtable.get(valueStr);
		if (cellBgColor != null)
		{
			cell.setBackground(cellBgColor);
		}
		Color cellFgColor = (Color) _cellFgColorHashtable.get(rowColVector);
		if (cellFgColor == null)
			cellFgColor = (Color) _cellFgColorHashtable.get(valueStr);
		if (cellFgColor != null)
		{
			cell.setForeground(cellFgColor);
        }

		// If this cell is selected, then highlight with a light-blue color
		if (isSelected) // && !(cell.getBackground().equals(Color.white)))
		{
			//cell.setBackground(cell.getBackground().darker());
			float[] rgb = cell.getBackground().getRGBComponents(null);
			rgb[0] *= .8;
			rgb[1] *= .8;  // darken the R&G components only, to highlight the blue component
			cell.setBackground(new Color(rgb[0],rgb[1],rgb[2]));
            rgb = cell.getForeground().getRGBComponents(null);
			rgb[0] *= .8;
			rgb[1] *= .8;  // darken the R&G components only, to highlight the blue component
            cell.setForeground(new Color(rgb[0],rgb[1],rgb[2]));
		}

		// If this cell should have a specific tooltip, then use it
		String cellTooltip = (String) _cellTooltipHashtable.get(rowColVector);
		if (cellTooltip == null)
			cellTooltip = (String) _cellTooltipHashtable.get(valueStr);
		if ((cellTooltip == null) || (cellTooltip.length() == 0))
		{
			// No specific tooltip set, so use the cell's string value as the tooltip
			if (value == null)
			{
				cell.setToolTipText(null);
			}
			else if (valueStr.length() > 200)
			{
				// Split long tooltip text into several smaller lines
				String tipText = "<html>";
				int MAX_CHARS_PER_LINE = 150;
				int strLen = valueStr.length();
				for (int lineIdx=0; lineIdx <= strLen/MAX_CHARS_PER_LINE; lineIdx++)
				{
					tipText = tipText.concat(valueStr.substring(lineIdx*MAX_CHARS_PER_LINE,Math.min((lineIdx+1)*MAX_CHARS_PER_LINE,strLen))).concat("<br>");
				}
				cell.setToolTipText(tipText);
			}
			else
			{
				cell.setToolTipText(valueStr);
			}
		}
		else
		{
			cell.setToolTipText(cellTooltip);
		}

		return cell;
	}

	// Constructors
	public ColoredFieldCellRenderer()
	{
		this(new Color(0.925F, 0.914F, 0.847F));  // gray
	}

	public ColoredFieldCellRenderer(Color bgcolor)
	{
		super();
		_bgcolor = bgcolor;
		//setOpaque(false);
		setBackground(_bgcolor);
	}

	public ColoredFieldCellRenderer(Color bgcolor, Color fgcolor)
	{
		super();
		_bgcolor = bgcolor;
		_fgcolor = fgcolor;
		//setOpaque(false);
		setBackground(_bgcolor);
		setForeground(_fgcolor);
	}

	public ColoredFieldCellRenderer(float[] bgrgb)
	{
		this(new Color(bgrgb[0], bgrgb[1], bgrgb[2]));
	}

	public ColoredFieldCellRenderer(float[] bgrgb, float[] fgrgb)
	{
		this(new Color(bgrgb[0], bgrgb[1], bgrgb[2]), 
             new Color(fgrgb[0], fgrgb[1], fgrgb[2]));
	}

	public ColoredFieldCellRenderer(float r, float g, float b)
	{
		super();
		_bgcolor = new Color(r,g,b);
		setBackground(_bgcolor);
	}

	// Default BG color
	public Color getBgColor()
	{
		return _bgcolor;
	}

	public void setBgColor(Color color)
	{
		_bgcolor = color;
		setBackground(_bgcolor);
	}

	public void setBgColor(float[] rgb)
	{
		_bgcolor = new Color(rgb[0], rgb[1], rgb[2]);
		setBackground(_bgcolor);
	}

	public void setBgColor(float r, float g, float b)
	{
		_bgcolor = new Color(r,g,b);
		setBackground(_bgcolor);
	}

	// Default FG color
	public Color getFgColor()
	{
		return _fgcolor;
	}

	public void setFgColor(Color color)
	{
		_fgcolor = color;
		setForeground(_fgcolor);
	}

	public void setFgColor(float[] rgb)
	{
		_fgcolor = new Color(rgb[0], rgb[1], rgb[2]);
		setForeground(_fgcolor);
	}

	public void setFgColor(float r, float g, float b)
	{
		_fgcolor = new Color(r,g,b);
		setForeground(_fgcolor);
	}

	// Reset BG colors / FG colors / Tooltips
	public void resetBgColors()
	{
        _cellBgColorHashtable.clear();
	}

	public void resetFgColors()
	{
        _cellFgColorHashtable.clear();
	}

	public void resetTooltips()
	{
        _cellTooltipHashtable.clear();
	}

	// Cell-specific BG Color / FG color / Tooltip
	public void setCellBgColor(int row, int column, Color color)
	{
		Vector rowColVector = getRowColVector(row, column);
		if (color == null)
		{
			color = Color.white;                         // Hashtables cannot accept nulls...
		}
		_cellBgColorHashtable.put(rowColVector, color);
	}

	public void setCellFgColor(int row, int column, Color color)
	{
		Vector rowColVector = getRowColVector(row, column);
		if (color == null)
		{
			color = Color.black;                         // Hashtables cannot accept nulls...
		}
		_cellFgColorHashtable.put(rowColVector, color);
	}

	public void setCellTooltip(int row, int column, String text)
	{
		Vector rowColVector = getRowColVector(row, column);
		if (text == null)
		{
			text = "";                                  // Hashtables cannot accept nulls...
		}
		_cellTooltipHashtable.put(rowColVector, text);
	}

	public Color getCellBgColor(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		return (Color) _cellBgColorHashtable.get(rowColVector);
	}

	public Color getCellFgColor(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		return (Color) _cellFgColorHashtable.get(rowColVector);
	}

	public String getCellTooltip(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		return (String) _cellTooltipHashtable.get(rowColVector);
	}

	private Vector getRowColVector(int row, int column)
	{
		Vector rowColVector = new Vector();
		rowColVector.addElement(new Integer(row));
		rowColVector.addElement(new Integer(column));
		return rowColVector;
	}

	// Value-specific BG Color / FG color / Tooltip
	public void setValueBgColor(String valueStr, Color color)
	{
		if (color == null)
		{
			color = Color.white;                         // Hashtables cannot accept nulls...
		}
		_cellBgColorHashtable.put(valueStr, color);
	}

	public void setValueFgColor(String valueStr, Color color)
	{
		if (color == null)
		{
			color = Color.black;                         // Hashtables cannot accept nulls...
		}
		_cellFgColorHashtable.put(valueStr, color);
	}

	public void setValueTooltip(String valueStr, String text)
	{
		if (text == null)
		{
			text = "";                                   // Hashtables cannot accept nulls...
		}
		_cellTooltipHashtable.put(valueStr, text);
	}

	public Color getValueBgColor(String valueStr)
	{
		return (Color) _cellBgColorHashtable.get(valueStr);
	}

	public Color getValueFgColor(String valueStr)
	{
		return (Color) _cellFgColorHashtable.get(valueStr);
	}

	public String getValueTooltip(String valueStr)
	{
		return (String) _cellTooltipHashtable.get(valueStr);
	}

	// Utility methods:
	// Debug
	public void setDebug(boolean flag)
	{
		_debug = flag;
	}

	public boolean isDebug()
	{
		return _debug;
	}

	// Disabled
	public void setDisabled(boolean flag)
	{
		_disabled = flag;
	}

	public boolean isDisabled()
	{
		return _disabled;
	}

	// SmartAlign
	public void setSmartAlign(boolean flag)
	{
		_smartAlign = flag;
	}

	public boolean isSmartAlign()
	{
		return _smartAlign;
	}
}
