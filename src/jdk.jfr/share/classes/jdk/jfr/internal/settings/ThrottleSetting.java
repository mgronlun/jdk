/*
 * Copyright (c) 2020, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.  Oracle designates this
 * particular file as subject to the "Classpath" exception as provided
 * by Oracle in the LICENSE file that accompanied this code.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 */

package jdk.jfr.internal.settings;

import java.util.Objects;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import jdk.jfr.Description;
import jdk.jfr.Label;
import jdk.jfr.MetadataDefinition;
import jdk.jfr.Name;
import jdk.jfr.Timespan;
import jdk.jfr.internal.PlatformEventType;
import jdk.jfr.internal.Type;
import jdk.jfr.internal.Utils;

@MetadataDefinition
@Label("Event Emmission Throttle")
@Description("Event emissions throttle (events per second)")
@Name(Type.SETTINGS_PREFIX + "Throttle")
public final class ThrottleSetting extends JDKSettingControl {
    private final static long typeId = Type.getTypeId(ThrottleSetting.class);
    private static final String OFF = "off";
    private static String DISABLED = "-1, -1.0, -1";

    private int rate = -1;
    private double probability = -1.0;
    private int nthSelection = -1;
    private String value = "-2, -1.0, -1";
    private final PlatformEventType eventType;

    public ThrottleSetting(PlatformEventType eventType) {
       this.eventType = Objects.requireNonNull(eventType);
    }

    @Override
    public String combine(Set<String> values) {
        long max = 0;
        String text = "0";
        System.out.println("Combine called");
        for (String value : values) {
            System.out.println("Combine value: " + value);
        }
        /*
        for (String value : values) {
            long l = parseValueSafe(value);
            if (l == 0) {
                // throttling is disabled, accept everything
                return "0";
            }
            if (l > max) {
                text = value;
                max = l;
            }
        }
        */
        return text;
    }

    @Override
    public void setValue(String value) {
        System.out.println("SetValue: " + value);
        this.value =  parseValueSafe(value);
        eventType.setThrottle(rate, probability, nthSelection);
    }

    @Override
    public String getValue() {
        return value;
    }

    public static boolean isType(long typeId) {
        return ThrottleSetting.typeId == typeId;
    }

    public String parseValueSafe(String value) {
        try {
            classifyAndParse(value);
        } catch (NumberFormatException nfe) {
        }
        if (rate == -1 && probability == -1.0 && nthSelection == -1) {
            // throttling is implicitly disabled; make it explicit by setting it to "off"
            rate = -2;
        }
        return rate + ", " + probability + ", " + nthSelection;
    }

    private void classifyAndParse(String value) {
        if (isDisabled(value)) {
            System.out.println("Throttler is " + value);
            rate = -2;
            return;
        }
        if (isRate(value)) {
            System.out.println(value + " is a rate");
            rate = (int)Utils.parseThrottleValue(value);
            return;
        }
        if (isProbability(value)) {
            System.out.println(value + " is a probability");
            parseProbability(value);
            return;
        }
        if (isNthSelection(value)) {
            System.out.println(value + " is an nth selection");
            parseNthSelection(value);
            return;
        }
        throw new NumberFormatException("'" + value + "' is not a valid input.");
    }

    private static boolean isDisabled(String value) {
        return OFF.equals(value.trim().toLowerCase());
    }

    private static boolean isRate(String value) {
        // Expected input format is "x / y" or "x \ y" where x is a non-negative integer and
        // y is a time unit.
        String trimmed = value.trim();
        if (trimmed.contains("/") && !trimmed.startsWith("/") && !trimmed.endsWith("/")) {
            return true;
        }
        return trimmed.contains("\\") && !trimmed.startsWith("\\") && !trimmed.endsWith("\\");
    }

    private static boolean isProbability(String value) {
        // Expected input format is "x %" or "0.x" where x is a non-negative integer or
        // a non-negative double value in the interval [0-1].
        return value.trim().endsWith("%") || value.contains(".");
    }

    private static boolean isNthSelection(String value) {
        // Expected input format is "x %" or "0.x" where x is a non-negative integer and
        // y is a time unit. Split the string at the delimiter.
        return value.trim().endsWith("th");
    }

    private void parseProbability(String value) {
        if (value.trim().endsWith("%")) {
            double parsedValue = Double.parseDouble(value.substring(0, value.length() - 1).trim());
            if (parsedValue < 0.0) {
                throw new NumberFormatException("'" + value + "' a probability can't be negative.");
            }
            if (parsedValue > 100) {
                probability = 1.0;
                return;
            }
            value = Double.toString(parsedValue / 100);
        }
        if (value.contains(".")) {
            String[] result = value.split("[.]");
            if (result.length != 2) {
                throw new NumberFormatException("'" + value + "' is not a valid format for a probability.");
            }
            result[0] = result[0].trim();
            result[1] = result[1].trim();
            System.out.println("Length: " + result.length);
            if (result[0] == "") {
                // empty characteristic, only mantissa ".y"
                System.out.println("Empty characteristic");
                probability = Double.parseDouble("0." + result[1]);
                return;
            }
            double characteristic = Double.parseDouble(result[0]);
            if (characteristic < 0) {
                throw new NumberFormatException("'" + value + "' a probability can't be negative.");
            }
            if (characteristic > 1.0) {
                probability = 1.0;
                return;
            }
            if (result[1] == "") {
                // empty mantissa, only characteristic "y."
                System.out.println("Empty mantissa");
                probability = Double.parseDouble(result[0] + ".0");
                return;
            }
            probability = Double.parseDouble(value.trim());
            return;
        }
        throw new NumberFormatException("'" + value + "' is not a valid format for a probability.");
    }

    private void parseNthSelection(String value) {
        nthSelection = Integer.parseInt(value.substring(0, value.length() - 2).trim());
    }
}
