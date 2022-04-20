const defaultTheme = require('tailwindcss/defaultTheme')
const colors = require('tailwindcss/colors')

module.exports = {
  content: [
    './layouts/**/*.html',
  ],
  theme: {
    extend: {
      typography: {
        DEFAULT: {
          css: {
            a: {
              color: '#104EB2',
            }
          },
        },
      },
      fontFamily: {
        sans: ['Open Sans', ...defaultTheme.fontFamily.sans]
      },
      colors: {
        'brand': {
          50: '#F3F6FB',
          100: '#E7EDF7',
          200: '#C3D3EC',
          300: '#9FB8E0',
          400: '#5883C9',
          500: '#104EB2',
          600: '#0E46A0',
          700: '#0A2F6B',
          800: '#072350',
          900: '#051735',
        },
        green: colors.emerald,
        yellow: colors.amber,
        purple: colors.violet
      }
    },
  },
  plugins: [
    require('@tailwindcss/typography'),
    require('@tailwindcss/line-clamp'),
  ],
}
