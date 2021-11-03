import { createMuiTheme } from "@material-ui/core/styles";
import publicSans from "../fonts/publicSans";

const theme = createMuiTheme({
  typography: {
    fontFamily: "Public Sans, sans-serif",
    h1: {
      fontFamily: "Merriweather",
      fontStyle: "normal",
      fontWeight: "bold",
      fontSize: "40px",
      color: "#000",
    },
    h2: {
      fontSize: "28px",
      color: "#000",
    },
    h3: {
      fontFamily: "Public Sans",
      fontStyle: "normal",
      fontWeight: "normal",
      fontSize: "22px",
      color: "#302E2E",
    },
    h4: {
      fontSize: "16px",
      fontWeight: "bold",
      color: "#000",
    },
    h5: {
      fontFamily: "Public Sans",
      fontStyle: "normal",
      fontWeight: "bold",
      fontSize: "16px",
      lineHeight: "24px",
      alignItems: "center",
      letterSpacing: "0.15px",
      color: "#000",
      display: "inline",
    },
    h6: {
      fontFamily: "Public Sans",
      fontStyle: "normal",
      fontSize: "14px",
      lineHeight: "24px",
      color: "#000",
      display: "inline",
    },
    subtitle1: {
      fontFamily: "Public Sans",
      fontStyle: "normal",
      fontWeight: "normal",
      fontSize: "14px",
      color: "#5C5C5C",
    },
    subtitle2: {
      fontFamily: "Public Sans",
      fontStyle: "normal",
      fontWeight: "normal",
      fontSize: "16px",
      lineHeight: "27px",
      display: "inline",
      color: "#000",
    },
  },
  palette: {
    primary: { main: "#005982", dark: "#162E51" },
    secondary: { main: "#007600", dark: "#006200" },
  },
  props: {
    MuiButtonBase: {
      disableRipple: true,
    },
  },
  overrides: {
    MuiCssBaseline: {
      "@global": {
        "@font-face": [publicSans],
      },
    },
    MuiTable: {
      stickyHeader: {
        borderCollapse: "collapse",
      },
    },
    MuiTableRow: {
      root: {
        padding: "0px",
      },
    },
    MuiTableCell: {
      root: {
        padding: "0px 16px",
        borderBottom: "0px",
      },
      stickyHeader: {
        zIndex: "auto",
      },
    },
    MuiTooltip: {
      tooltip: {
        backgroundColor: "#747474",
        borderRadius: "8px",
        color: "#FFFFFF",
        fontSize: "16px",
      },
    },
    MuiCheckbox: {
      colorPrimary: {
        color: "#535353",
        "&$checked": {
          color: "#005EA2",
        },
        "&$disabled": {
          color: "#979797",
        },
      },
    },
    MuiRadio: {
      colorPrimary: {
        color: "#535353",
        "&$checked": {
          color: "#005EA2",
        },
        "&$disabled": {
          color: "#979797",
        },
      },
    },
    MuiChip: {
      colorPrimary: {
        color: "#FFFFFF",
        backgroundColor: "#005982",
      },
      outlinedPrimary: {
        color: "#535353",
        border: "1px solid #535353",
      },
    },
    MuiButton: {
      label: {
        fontWeight: "bold",
        textTransform: "none",
      },
      containedPrimary: {
        backgroundColor: "#06ac0c",
        "&:hover": {
          boxShadow: "0 0 0 0",
          backgroundColor: "#007600",
        },
        "&$disabled": {
          color: "white",
          backgroundColor: "#949494",
        },
        "&:active, &$focusVisible": {
          boxShadow: "0 0 0 0",
          backgroundColor: "#007600",
        },
      },
      containedSecondary: {
        boxShadow: "0 0 0 0",
        "&:active": {
          boxShadow: "0 0 0 0",
          backgroundColor: "#006200",
        },
        "&:hover": {
          boxShadow: "0 0 0 0",
          backgroundColor: "#006200",
        },
      },
      outlinedSecondary: {
        boxShadow: "0 0 0 0",
        color: "#007600",
        borderWidth: "2px !important",
        "&:active, &$focusVisible": {
          boxShadow: "0 0 0 0",
          color: "white",
          backgroundColor: "#006200",
        },
        "&:hover": {
          boxShadow: "0 0 0 0",
          color: "white",
          backgroundColor: "#007600",
        },
      },
    },
    MuiSvgIcon: {
      // secondary = success = green
      colorSecondary: {
        color: "#007600",
      },
      // action = info = orange
      colorAction: {
        color: "#F4AD22",
      },
      // error = warning = red
      colorError: {
        color: "#C20B0B",
      },
      fontSizeSmall: {
        fontSize: "1rem",
      },
      fontSizeLarge: {
        fontSize: "1.5rem",
      },
    },
    MuiLink: {
      root: {
        fontWeight: "bold",
      },
    },
    MuiAccordion: {
      root: {
        "&$expanded": {
          margin: 0,
          borderRight: "1px solid #00829B",
          borderLeft: "1px solid #00829B",
        },
      },
    },
    MuiAlert: {
      icon: {
        padding: 0,
      },
      message: {
        padding: 0,
      },
      filledSuccess: {
        backgroundColor: "#D9FCDA",
        color: "#949494",
      },
      filledWarning: {
        backgroundColor: "#FFE9BD",
        color: "#949494",
      },
      filledInfo: {
        backgroundColor: "#BFF4FF",
        color: "#949494",
      },
      filledError: {
        backgroundColor: "#FFE4DE",
        color: "#949494",
      },
    },
    MuiAvatar: {
      root: {
        fontSize: "1.0rem",
        width: "30px",
        height: "30px",
      },
    },
    MuiList: {
      root: {
        outline: "none",
      },
    },
  },
});

export default theme;
